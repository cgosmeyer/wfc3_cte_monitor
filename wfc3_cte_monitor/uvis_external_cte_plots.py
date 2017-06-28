"""
Module containing functions that create the primary plots for the 
analysis of the WFC3/UVIS External CTE monitor.

Authors:

    C.M. Gosmeyer, Spring 2016

Use:

    Nominally functions are imported to ``run_outputs.py``. Can also be
    run from the command line. You will need adjust function calls and
    parameters in the Main.

    >>> python uvis_external_cte_plots.py
 
"""

import itertools
import numpy as np
import os
import pylab

from collections import OrderedDict
from scipy.stats import linregress
from scipy.stats import sigmaclip

from config import program_list
from config import filter_list
from config import target_list
from config import flashlvl_list
from config import fluxbins_lo
from config import fluxbins_hi

from analysis_tools.statistics.meanclip import meanclip

from database_queries import return_tables
from database_queries import query_for_pair
from database_queries import query_for_flux_by_imagename
from database_queries import query_for_exptimes
from database_queries import query_for_dateobss
from database_queries import query_results_for_slopes
from database_queries import query_for_fluxes_bkgrds_by_ypos
from database_queries import query_for_globalbkgrd
from database_queries import query_for_all_dateobss
from database_queries import query_for_flux_range
from database_queries import query_results_for_slopes
from database_queries import query_for_matching_imagename
from database_queries import query_for_180pair

from database_update import return_results_dict
from database_update import update_results_table

from database_interface import NGC104_Results
from database_interface import NGC6583_Results
from database_interface import NGC6791_Results


#-------------------------------------------------------------------------------# 
# Base plot functions.
#-------------------------------------------------------------------------------# 

def plot_fluxratio_vs_ypos(imagename_chip1, imagename_chip2, fluxes_chip1, 
                     fluxes_chip2, bkgrds_chip1, bkgrds_chip2, ypos_chip2,
                     master_id_chip2, param_dict, chip=0, 
                     write_ratios_to_file=False, outloc=''):
    """ 
    Ratio of the fluxes (chip1/chip2) over chip2's y-position for a 
    single image pair over various flux bins. The procedure is
        * Subtract annulus backgrounds from flux.
        * Sigma clip the flux ratio once.
        * Take the linear fit and returns a slope. 
        * Assume that have already weeded out the non-matches of flux 1
            and 2 using the master_id. 
        * If not enough points in a flux bin to plot, then insert NaNs 
            into the databse for the slope and slope standard dev of that 
            flux bin. 
            (Perhaps NaNs are not the accepted standard for database
            entries. you might consider changing this.)
    Creates a seperate plot and thus slope for each flux bin. Each flux 
    bin's slope gets a column in the database for the given image pair's 
    row.

    If the 180-degree dataset (specified if set 'chip=0') then the ratio 
    is first visit / second visit vs y-position of second visit.

    Parameters:
        imagename_chip1 : string
            Name of the image containing chip1. (Or first visit if 180-
            degree dataset.)
        imagename_chip2 : string
            Name of the image containing chip2. (Or second visit if 180-
            degree dataset.)
        fluxes_chip1 : array
            Fluxes of chip1 sources.
        fluxes_chip2 : array
            Fluxes of chip2 sources.
        bkgrds_chip1 : array
            Local backgrounds of chip1 sources.
        bkgrds_chip2 : 
            Local backgrounds of chip2 sources.
        ypos_chip2 : array
            Y-positions of chip2 sources.
        master_id_chip2 : array
            Master IDs of chip2 sources.
        param_dict : dict
            Dictionary containing image's header values and other
            useful information that will be used to fill in database
            and label plots and data files.
        chip : int
            If querying for 180-degree datasets, set this parameter for 
            the chip number of desired pair. Otherwise set to 0.
        write_ratios_to_file : {True, False}
            Set to True to record the ratios for each flux bin in a file 
            '*_slopes.txt'.
        outloc : string
            Path to the output location. Working directory by default.

    Returns:
        String. Name of plot.

    Outputs:
        PNG plot. If nominal dataset,
        "<imagename1>_<imagename2>_r<aperture>_slopes.png"
        And if 'write_ratios_to_file' set to True, text file,
        "<imagename1>_<imagename2>_r<aperture>_slopes.txt"

        If 180-degree dataset,
        "<imagename1>_<imagename2>_r<aperture>_ch<1/2>_slopes.png"    
        And if 'write_ratios_to_file' set to True, text file,
        "<imagename1>_<imagename2>_r<aperture>_ch<1/2>_slopes.txt"        

    Notes:
        * For the 180-degree data, simply put the first and second
          visit images in place of chip1 and chip2.

    References:
        http://stackoverflow.com/questions/6148207/linear-regression-with-matplotlib-numpy
    """
    # Subtract the annulus backgrounds.
    fluxes_chip1_clean = np.array(fluxes_chip1) - np.array(bkgrds_chip1) 
    fluxes_chip2_clean = np.array(fluxes_chip2) - np.array(bkgrds_chip2) 

    # Sort fluxes into bins.
    # Each shall be its own color. 
    # Loop through each bin, removing outliers, fitting a line, and 
    # writing slopes to a file.
    slopes = []
    slope_stddevs = []
    num_points_ls = []
    colors = ['grey', 'indigo', 'blue', 'green', 'brown', 'orange', 
              'red', 'magenta']
    labels = ['250-500', '500-1000', '500-2000', '1000-2000', 
              '2000-4000', '2000-8000', '4000-8000', '8000-32000']

    # Set plot.
    fig=pylab.figure(figsize=(12.5,8.5))

    # Loop over the flux bins.
    for k in range(len(labels)):
        # Weed out entries that don't fit in the bin.
        binned_indices = [x[0] for x,y in zip(enumerate(fluxes_chip1_clean), 
                                               enumerate(fluxes_chip2_clean)) 
                          if ((x[1] > fluxbins_lo[k] and x[1] < fluxbins_hi[k]) and 
                              (y[1] > fluxbins_lo[k] and y[1] < fluxbins_hi[k]))]     
        fluxes_chip1_binned = fluxes_chip1_clean[binned_indices] 
        fluxes_chip2_binned = fluxes_chip2_clean[binned_indices] 
        master_id_chip2_binned = master_id_chip2[binned_indices]
        ypos_chip2_binned = ypos_chip2[binned_indices]

        # Calculate the flux ratios of all remaining entries.
        flux_ratio = fluxes_chip1_binned/fluxes_chip2_binned

        # Sigma clip to remove outliers.
        flux_ratio_sigclipped, thresh_low, thresh_up = \
            sigmaclip(flux_ratio, high=4, low=4)

        # Match the indices of sig clipped array back to ypos array.
        matched_indices = match_indices(flux_ratio_sigclipped, flux_ratio) 
        master_id_chip2_sigclipped = master_id_chip2_binned[matched_indices]
        ypos_chip2_sigclipped = ypos_chip2_binned[matched_indices]
        if write_ratios_to_file:
            bkgrds_chip1_sigclipped = bkgrds_chip1[matched_indices]
            bkgrds_chip2_sigclipped = bkgrds_chip2[matched_indices]
            fluxes_chip1_sigclipped = fluxes_chip1[matched_indices]
            fluxes_chip2_sigclipped = fluxes_chip2[matched_indices]
        
        # Fit line to binned, sigma-clipped fluxes. 
        if ypos_chip2_sigclipped != [] and flux_ratio_sigclipped != []:
            m, b, r_value, p_value, std_err = linregress(ypos_chip2_sigclipped, 
                flux_ratio_sigclipped)
            line = m*np.array([np.min(ypos_chip2), np.max(ypos_chip2)]) + b
            # Calculate standard deviation from standard error, where
            # stddev = stderr * srt(number of points)
            num_points = len(ypos_chip2_sigclipped)
            std_dev = std_err * np.sqrt(num_points)

        else:
            m = 0
            num_points = 0
            std_dev = 0

        slopes.append(m)
        slope_stddevs.append(std_dev)
        num_points_ls.append(num_points)

        if write_ratios_to_file:
            # create file name
            if chip == 0:
                ratios_filename = '{}_{}_{}_r{}_fluxratios.txt'\
                        .format(imagename_chip1.split('.fits')[0], 
                            imagename_chip2.split('.fits')[0], 
                            labels[k], param_dict['aperture'])
            else: # for 180-degree dataset, label chip used
                 ratios_filename = '{}_{}_{}_r{}_ch{}_fluxratios.txt'\
                        .format(imagename_chip1.split('.fits')[0], 
                            imagename_chip2.split('.fits')[0], 
                            labels[k], param_dict['aperture'], chip)               
            print("Writing ratios to a *fluxratios.txt file...")
            # Write the flux ratios to a file.
            # Check first whether the file already exists. If so, append.
            if os.path.isfile(os.path.join(outloc, ratios_filename)):
                open_file = open(os.path.join(outloc, ratios_filename), 'a')
            # If it doesn't exist, create a header. 
            else: 
                open_file = open(os.path.join(outloc, ratios_filename), 'w')
                open_file.write("#{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                    .format('master_id', 'fluxratio_sigclpd', 
                    'flux_chip1', 'flux_chip2', 'bkgrd_chip1', 
                    'bkgrd_chip2', 'ypos') )

            for masterid, fluxratio, flux_chip1, flux_chip2, bkgrd_chip1, bkgrd_chip2, ypos in \
                zip(master_id_chip2_sigclipped, flux_ratio_sigclipped, 
                    fluxes_chip1_sigclipped, fluxes_chip2_sigclipped, 
                    bkgrds_chip1_sigclipped, bkgrds_chip2_sigclipped, ypos_chip2_sigclipped):
                open_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(masterid, fluxratio, 
                    flux_chip1, flux_chip2, bkgrd_chip1, bkgrd_chip2, ypos) )
            open_file.close()

        # Only plot if there are more than two entries in the bin. Otherwise
        # line fitter explodes.
        if num_points >= 3:
            # Draw lines at ratio=1.0 and y-postion=1024
            pylab.axhline(y=1.0, linestyle='--', color='black')
            pylab.axvline(x=1024, linestyle='--', color='black')
            # Plot the flux ratio vs ypos.
            pylab.scatter(ypos_chip2_sigclipped, flux_ratio_sigclipped, 
                s=10, marker='x', color=colors[k], label=labels[k]) 
            # draw the slope line
            pylab.plot([np.min(ypos_chip2_sigclipped), np.max(ypos_chip2_sigclipped)], 
                line, '-', color=colors[k], linewidth='1.5')  
        else:
            print("Skipping the plotting of bin {}-{} because only 2 or fewer points."\
                .format(fluxbins_lo[k], fluxbins_hi[k]))
            print("Inserting NaNs into database for slope and slopestdev.")

    # Complete the plot and save. 
    if chip == 0:
        pylab.xlabel('Y-position_chip2 [pxl]', fontsize=22, weight='bold') 
        pylab.ylabel('Flux_chip1 / Flux_chip2', fontsize=22, weight='bold')
        title = '{} {} mjd{} {} {}s pf{} cte{} ap{}'\
           .format(param_dict['targname'], param_dict['proposid'],
            param_dict['dateobs'], param_dict['filter'], 
            param_dict['exptime'], param_dict['flashlvl'], 
            param_dict['ctecorr'], param_dict['aperture'])
    else: # for 180-degree dataset, label chip used
        pylab.xlabel('Y-position_visit2 [pxl]', fontsize=22, weight='bold') 
        pylab.ylabel('Flux_visit1 / Flux_visit2', fontsize=22, weight='bold')
        title = '{} {} mjd{} {} {}s pf{} cte{} ap{} ch{}'\
           .format(param_dict['targname'], param_dict['proposid'],
            param_dict['dateobs'], param_dict['filter'], 
            param_dict['exptime'], param_dict['flashlvl'], 
            param_dict['ctecorr'], param_dict['aperture'], chip)
    pylab.tick_params(axis='both', which='major', labelsize=20)
    pylab.title(title, fontsize=16)
    pylab.xlim([0,2000])
    pylab.ylim([0.2,1.6])
    pylab.legend(scatterpoints=3, loc='best')

    # Create outname for plot and text file.
    if chip == 0:
        outname = '{}_{}_r{}'.format(imagename_chip1.split('.fits')[0], 
                               imagename_chip2.split('.fits')[0], 
                               param_dict['aperture'])
    else: # for 180-degree dataset, label chip used
        outname = '{}_{}_r{}_ch{}'.format(imagename_chip1.split('.fits')[0], 
                               imagename_chip2.split('.fits')[0], 
                               param_dict['aperture'], chip)        
    pylab.savefig(os.path.join(outloc, '{}_slopes.png'.format(outname)), 
        bbox_inches='tight')


    # Write the slopes and related info to text file.
    # First write in a header with param_dict.
    open_file = open(os.path.join(outloc, '{}_slopes.txt'.format(outname)), 'w')
    for key, val in zip(param_dict.keys(), param_dict.values()):
        open_file.write('#{}\t{}\n'.format(key,val))
    open_file.write
    open_file.write("#{}\t{}\t{}\t{}\t{}\n".format('slope', 'slope_stddev', 'num_points',
        'low_bin', 'high_bin') )
    for slope, slope_stddev, num_points, low_bin, high_bin in zip(slopes, 
        slope_stddevs, num_points_ls, fluxbins_lo, fluxbins_hi):
        open_file.write("{}\t{}\t{}\t{}\t{}\n".format(slope, slope_stddev, num_points,
            low_bin, high_bin) )
    open_file.close()

    return imagename_chip1, imagename_chip2, '{}_slopes.txt'.format(outname)


#-------------------------------------------------------------------------------# 

def plot_cteslope_vs_time(slopes, times, standerrs, param_dict, outloc='', 
    xlims=[55000,57800], ylims=[-0.1, 0.6], write_to_file=True):
    """ Chip1/Chip2 flux ratio slopes (CTE measure) over time for a postflash-
    exptime-filter-ctecorr-aperture group.  Overplots ngc104 and ngc6791 
    as different marker types and three different flux bins as different 
    colors.

    This is the primary plot for analyzing CTE evolution with time.

    Parameters:
        slopes : list of arrays
            The flux-ratio vs y-position slopes of a parameter set.
        times : list of arrays
            The dates for each slope.
        standerrs : list of arrays
            The standard error for each slope.
        param_dict : dict
            Dictionary containing image's header values and other
            useful information that will be used to fill in database
            and label plots and data files.
        outloc : string
            Path to the output location. Working directory by default.
        xlims : list of floats
            The lower and upper bounds of the x-axis (MJD time). 
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope). 
        write_to_file : {True, False}
            Set on to write the slope and time lists to a text file.
            True by default.

    Returns:
        String. Name of the plot.

    Outputs:
        PNG plot. 
        "cteVStime_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>.png"
    """
    # Create outname for plot and text file.
    outname = 'cteVStime_{}_{}_pf{}_ctecorr{}_r{}'.format(param_dict['filter'], 
        param_dict['exp_length'], param_dict['flashlvl'], 
        param_dict['ctecorr'], param_dict['aperture'])

    # Match the markers and colors with appropriate targnames and filter bins.
    markers = {'ngc104' : 'o', 'ngc6791' : 'x'}
    colors = {'500-2000' : 'k', '2000-8000' : 'red', '8000-32000' : 'blue'}

    # Set the figure.
    pylab.figure(figsize=(12.5,8.5))

    # Loop over targets.
    for targname in slopes.keys():
        # Loop over fluxbins. 
        for fluxbin in (slopes[targname]).keys():

            # Change all Nones to NaNs.
            slopes_bin = np.array((slopes[targname])[fluxbin], dtype=np.float)
            if standerrs != []:
                standerrs_bin = np.array((standerrs[targname])[fluxbin], dtype=np.float)
            times_bin = np.array((times[targname])[fluxbin], dtype=np.float)

            # Correct the CTE measurement for the chip-height (2048 pixels)
            slopes_bin_corrected = (slopes_bin / 2.) * 2048.
            standerrs_bin_corrected = (standerrs_bin / 2.) * 2048.

            # Draw the scatter plot for target/fluxbin combination.
            pylab.scatter(times_bin, slopes_bin_corrected, 
                s=90, marker=markers[targname], color=colors[fluxbin], 
                label='{} {}e-'.format(targname, fluxbin), alpha=0.6)
            if standerrs != []:
                pylab.errorbar(times_bin, slopes_bin_corrected,  
                    yerr=standerrs_bin_corrected, 
                    color='gray', alpha=0.75, marker=None, ls='None')

            # Write slope and time lists to text file.
            if write_to_file:
                file_name = os.path.join(outloc, '{}.txt'.format(outname))
                if os.path.isfile(file_name):
                    # If already exists, just append.
                    open_file = open(file_name, 'a')
                else:
                    # If doesn't exist already, create header.
                    open_file = open(file_name, 'w')
                    open_file.write("#{}\t{}\t{}\t{}\t{}\n".format('slope', 
                        'slope_stderr', 'mjd', 'targname', 'fluxbin') )
                
                # Also add targname and fluxbin into file.
                for slope, slope_stderr, mjd in zip(slopes_bin_corrected, 
                    standerrs_bin_corrected, times_bin):
                    open_file.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\n".format(slope, 
                        slope_stderr, mjd, targname, fluxbin) )
                open_file.close()

    # Complete the plot and save. 
    pylab.xlabel('MJD Date', fontsize=22, weight='bold') 
    pylab.ylabel('CTE loss [flux / 2048 pxl]', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    # Add lines at major y tickmarks. Must define by hand...
    for ymaj in np.arange(11)*0.1:
        pylab.axhline(y=ymaj,linestyle='--',color='grey')
    title = "{} explen'{}' pf{} ctecorr{} ap{}".format(param_dict['filter'], 
        param_dict['exp_length'], param_dict['flashlvl'], 
        param_dict['ctecorr'], param_dict['aperture'])
    pylab.title(title, fontsize=16)
    pylab.xlim(xlims)
    pylab.ylim(ylims)
    pylab.legend(scatterpoints=1, loc='upper left')
    if standerrs != []:
        pylab.text(56500, 0.55, 'Standard Error: stddev / sqrt(# sources)')

    pylab.savefig(os.path.join(outloc, '{}.png'.format(outname)), 
        bbox_inches='tight')

    return '{}.png'.format(outname)


#-------------------------------------------------------------------------------# 

def plot_cteslope_vs_logflux(targname, slopes, flux_means, all_fluxes, 
    standerrs, param_dict, chip=0, outloc='', ylims=[-0.1,0.6], 
    skipbins=0, clipoutliers=False, 
    plot_errorbars=False, chip1_invert=True):
    """ Chip1/Chip2 flux ratio slopes (CTE measure) vs the log10flux 
    of the average of the bin (e.g., the 500-1000 bin will average to 750). 
    Uses measurements from photometry aperture 3 radius. Change this by
    retrieving the 'slopes', 'flux_means', etc from a different radius
    stored in the database. 
    Fits a 2nd degree polynomial to each epoch and returns the 9 
    coefficients to the fit.

    The primary plot for analyzing CTE evolution with flux strength. 

    Parameters:
        targname : string
            Name of the target cluster.
        slopes : dict of lists
            Keys of epoch. Values of lists of CTE slopes.
        flux_means : dict of lists
            Keys of epoch. Values of lists of mean fluxes.
        all_fluxes : dict of lists
            Keys of epoch. Values of lists of once-sigma-clipped fluxes.
        standerrs : 
            The standard error for each slope.
        param_dict : dict
            Dictionary containing image's header values and other
            useful information that will be used to fill in database
            and label plots and data files.
        chip : int
            If querying for 180-degree datasets, set this parameter for 
            the chip number of desired pair. Otherwise set to 0.
        outloc : string
            Path to the output location. Working directory by default.
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope). 
        skipbins : int
            Number of the bins to skip. Starts at 1. O indicates no skipping.
        clipoutliers : {True, False}
            Turn on to clip 2*stddev outliers from polynomial fit.
            Not very affective. Recomended to play with 'skipbins' instead.
        plot_errorbars : {True, False}
            Set on to plot errorbars.
        chip1_invert : {True, False}
            Multipy -1.0 to values of the chip1 fluxes. True by default.
            (180-degree dataset's chip1 slopes negative because direction of 
            readout in opposite direction from chip2's direction.)

    Returns:
        String. Name of the plot.

    Outputs:
        PNG plot. If nominal dataset,
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>.png"

        If 180-degree dataset,
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>_ch<1/2>.png"

        Text files containing polynomial fit's 9 coefficients for each
        epoch.
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>_mjd<epoch>_coeffs.txt"
    """

    # Create outname for plot and text file.
    if chip == 0:
        outname = 'cteVSlogflux_{}_{}_{}_pf{}_ctecorr{}_r{}'.format(
            targname, param_dict['filter'], param_dict['exp_length'], 
            param_dict['flashlvl'], param_dict['ctecorr'],
            param_dict['aperture'])
    else:
        outname = 'cteVSlogflux_{}_{}_{}_pf{}_ctecorr{}_r{}_ch{}'.format(
            targname, param_dict['filter'], param_dict['exp_length'], 
            param_dict['flashlvl'], param_dict['ctecorr'],
            param_dict['aperture'], chip)

    # One color for each epoch, which cycles using itertools.
    colors_ls = ['brown', 'red',  'gold', 'lime', 'green', 
                 'cyan', 'blue', 'purple', 'magenta', 'black']
    colors = itertools.cycle(colors_ls)

    # Set plot.
    fig=pylab.figure(figsize=(12.5,8.5))
    
    # Set a count, after which start plotting lines in different 
    # style, so that no two epochs have same color and style.
    epoch_count = 0

    # Loop over available epochs.
    for epoch in slopes.keys():
        # For each epoch and flux bin, need to average all slopes.
        # Take the log10 of the average flux of the flux bin.
        if flux_means[epoch] != []:
            if epoch_count < len(colors_ls):
                linestyle = '--'
            else:
                linestyle = '-'
            epoch_count += 1

            # Change all Nones to NaNs.
            if chip == 1 and chip1_invert:
                slopes_epoch = (-1.0)*(np.array(slopes[epoch], dtype=np.float))
            else:
                slopes_epoch = np.array(slopes[epoch], dtype=np.float)
            flux_means_epoch = np.array(flux_means[epoch], dtype=np.float)

            # Set the next color in sequence.
            color = next(colors)

            # Get x and y values into desired formats.
            log10fluxes = np.log10(flux_means_epoch)
            corrected_cteslopes = (slopes_epoch / 2.) * 2048.

            # Plot!
            pylab.scatter(log10fluxes, corrected_cteslopes, 
                s=90, marker='d', alpha=0.75, color = color)

            if plot_errorbars:
                standerrs_epoch = np.array(standerrs[epoch], dtype=np.float)
                corrected_standerrs = (standerrs_epoch / 2.) * 2048.
                pylab.errorbar(log10fluxes, corrected_cteslopes,
                    yerr=corrected_standerrs, 
                    color='gray', alpha=0.75, marker=None, ls='None')

            # Overplot the empirical model fit.
            model_fit = fit_empirical_model(log10fluxes=log10fluxes, 
                                            corrected_cteslopes=corrected_cteslopes, 
                                            epoch=epoch,
                                            outname=outname,
                                            outloc=outloc,
                                            all_fluxes_epoch=all_fluxes[epoch],
                                            skipbins=skipbins,
                                            clipoutliers=clipoutliers)
            pylab.plot(model_fit[0], model_fit[1], color=color, 
                linewidth=2.5, linestyle=linestyle, label='MJD {}'.format(int(epoch)))


    # Complete the plot and save. 
    pylab.xlabel('LOG Flux (ap 3 pxl radius) [e-]', fontsize=22,
        weight='bold') 
    pylab.ylabel('CTE loss [flux / 2048 pxl]', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    # Add lines at major y tickmarks. Must define by hand...
    for ymaj in np.arange(11)*0.1:
        pylab.axhline(y=ymaj,linestyle='--',color='grey')
    if chip == 0:
        title = "{} {} explen'{}' pf{} ctecorr{} ap{}".format(
            targname, param_dict['filter'], 
            param_dict['exp_length'], param_dict['flashlvl'], 
            param_dict['ctecorr'], param_dict['aperture'])
    else:
        title = "{} {} explen'{}' pf{} ctecorr{} ap{} ch{}".format(
            targname, param_dict['filter'], 
            param_dict['exp_length'], param_dict['flashlvl'], 
            param_dict['ctecorr'], param_dict['aperture'], chip)        
    pylab.title(title, fontsize=16)
    pylab.xlim([2.5, 4.5])
    pylab.ylim(ylims)
    pylab.legend(scatterpoints=1, loc='upper right', ncol=3)

    pylab.savefig(os.path.join(outloc, '{}.png'.format(outname)), 
        bbox_inches='tight')

    return '{}.png'.format(outname)


#-------------------------------------------------------------------------------# 

def fit_empirical_model(log10fluxes, corrected_cteslopes, epoch, outname, 
    outloc='', all_fluxes_epoch=[], skipbins=0,
    clipoutliers=False, write_fit_values=False, chip=0):
    """ Fits a 2-D polynomial to an epoch's cte-slope vs log flux.
    The function `polyfit2d` is defined in this module. The original IDL code 
    used the procedure `sfit`.

    S : CTE slope. [e-/2048 pxl]
    d : observation date, MJD - 55400
    f : source flux, log10(flux[e-])
    S = c_00 + c_01*f + c_02*f^2 + c_10*d + c_11*f*d + c_12*f^2*d + 
        c_20*d^2 + c_21*f*d^2 + c_22*f^2*d^2
        (note that this equation was incorrect in WFC3 ISRs 2012-09 and 
        2015-03, but was always correct in the scripts)

    Do a fit for each epoch. Users should be using coefficients from the
    latest epoch. 

    Parameters:
        log10fluxes : array of floats
            Fluxes (average of flux bins) at log10.
        corrected_cteslopes : array of floats
            CTE slopes, corrected for distance from amp.
        epoch : float or int
            The epoch of the slopes in MJD.
        outname : string
            Name of the output text file.
        outloc : string
            Path to the output location. Working directory by default.
        all_fluxes_epoch : dict of lists
            Keys of epoch. Values of lists of fluxes.
        skipbins : int
            Number of the bins to skip. Starts at 1. O indicates no skipping.
        clipoutliers : {True, False}
            Turn on to clip 2*stddev outliers from polynomial fit.
            Not very affective. Recomended to play with 'skipbins' instead.
        write_fit_values : {True, False}
            Turn on to write the fit values for the cteloss and logflux axis
            to a text file "*fitvals.txt".

    Returns:
        fluxgrid : array of floats
            Interpolated log10 fluxes from 2.5 to 4.48. (x-axis)
        fitvals : array of floats
            Values of 2-D polynomial fit over the 'fluxgrid'. (y-axis)

    Outputs:
        Text files containing polynomial fit's 9 coefficients for each
        epoch.
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>_mjd<epoch>_coeffs.txt"

        If 'write_fit_values' set to True,
        Text files containing the logflux vs time fit values.
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>_mjd<epoch>_fitvals.txt"        
    
    References:
        http://stackoverflow.com/questions/10988082/multivariate-polynomial-regression-with-numpy
        http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
    """
    # Define the poly degree. (the "2-d" bit)
    polydeg = 2

    # Fill out the model values
    log10fluxes = np.array(log10fluxes)
    corrected_cteslopes = np.array(corrected_cteslopes)
    epochs = np.ones(len(log10fluxes))*(epoch-55400.)

    # Because fit seems to be thrown off by low flux bins, provide option
    # to skip them.
    log10fluxes_0 = log10fluxes[skipbins:] 
    corrected_cteslopes_0 = corrected_cteslopes[skipbins:]
    epochs_0 = epochs[skipbins:]

    # NaNs are now breaking polyfit2d in Python 3 so get rid of them.
    epochs = []
    log10fluxes= []
    corrected_cteslopes = []
    for mjd, log10flux, corrected_cteslope in zip(epochs_0, log10fluxes_0, corrected_cteslopes_0):
        if not np.isnan(corrected_cteslope):
            print(corrected_cteslope)
            epochs.append(mjd)
            log10fluxes.append(log10flux)
            corrected_cteslopes.append(corrected_cteslope)
        print(corrected_cteslopes)

    epochs = np.array(epochs)
    log10fluxes = np.array(log10fluxes)
    corrected_cteslopes = np.array(corrected_cteslopes)

    print(epochs)
    print(log10fluxes)
    print(corrected_cteslopes)

    # Get 9 the coefficients.
    coeffs = polyfit2d(epochs, log10fluxes, corrected_cteslopes)
    print(coeffs)

    # Initialize a grid on which to evaluate coeffs.
    fluxgrid = (np.arange(100))/50. + 2.5

    # Evaluate them.
    fitvals = np.zeros(len(fluxgrid))
    coeff_count = 0
    for row in range(polydeg+1):
        for col in range(polydeg+1):
            print('coeff_count: {}'.format(coeff_count))
            fitvals += coeffs[coeff_count]*(fluxgrid**col)*((epoch-55400.)**row)
            print('fitvals: {}'.format(fitvals))
            coeff_count += 1

    # Perform second polyfit, removing outliers.
    if clipoutliers:
        # Actually this does not do much. SKipping bins is more affective. 
        # What really needs be done is more clipping in the CTE slopes themselves.
        print("Clipping outliers...")
        flux_means_clip = []
        for fluxes in all_fluxes_epoch:
            print('length fluxes:', len(fluxes))
            fluxes_clip = sigmaclip(fluxes, high=1.3, low=1.3)[0]
            print('length fluxes_clip:', len(fluxes_clip))
            flux_mean = np.mean(fluxes_clip)
            flux_means_clip.append(flux_mean)
        log10fluxes_clip = np.log10(flux_means_clip)
        
        # Skip desired number of bins
        log10fluxes_clip = log10fluxes_clip[skipbins:]

        print(len(log10fluxes_clip))
        print(len(epochs))
        print(len(corrected_cteslopes))
        print("orig log10fluxes: {}".format(log10fluxes))
        print("clipped log10fluxes: {}".format(log10fluxes_clip))

        # Get 9 coefficients
        coeffs = polyfit2d(epochs, log10fluxes_clip, corrected_cteslopes)
        print(coeffs)

        # Evaluate them on a grid again.
        fluxgrid = (np.arange(100))/50. + 2.5

        fitvals = np.zeros(len(fluxgrid))
        coeff_count = 0
        for row in range(polydeg+1):
            for col in range(polydeg+1):
                print('coeff_count: {}'.format(coeff_count))
                fitvals += coeffs[coeff_count]*(fluxgrid**col)*((epoch-55400.)**row)
                print('fitvals: {}'.format(fitvals))
                coeff_count += 1

    # Print to a file *coeffs.dat.
    # Create file name
    coeff_filename = '{}_mjd{}_coeffs.txt'.format(outname, int(epoch))
    print("Writing coeffs to a *coeffs.txt file...")
    # Write the coeffs to a file.
    open_file = open(os.path.join(outloc, coeff_filename), 'w')
    open_file.write("#polydeg: {}\n".format(polydeg))
    open_file.write("#used Python function 'polyfit2d'\n")
    open_file.write("#{}\t{}\t{}\n".format('row', 'col', 'coeff'))

    # Write the coefficients.
    i = 0
    for row in range(polydeg+1):
        for col in range(polydeg+1):
            open_file.write("{}\t{}\t{}\n".format(row, col, coeffs[i]) )
            i += 1
    open_file.close()

    # Write the fit values if option set to True.
    if write_fit_values:
        if chip == 0:
            # Nominal dataset.
            fit_filename = '{}_mjd{}_fitvals.txt'.format(outname, int(epoch))
        else:
            # 180degree dataset
            fit_filename = '{}_mjd{}_ch{}_fitvals.txt'.format(outname, int(epoch), chip)
        print("Writing logflux and time fit values to a *fitvals.txt file...")
        open_file = open(os.path.join(outloc, fit_filename), 'w')
        open_file.write("#{}\t{}\n".format('cteloss', 'logflux'))
        for cteloss, logflux in zip(fitvals, fluxgrid):
            open_file.write("{}\t{}\n".format(cteloss, logflux))
        open_file.close()

    return fluxgrid, fitvals


#-------------------------------------------------------------------------------# 

def plot_flux_vs_bkgrd(image_dict, param_dict, local_bkgrd=True, outloc=''):
    """ This should emulate Figure 3 of Jay's ISR. 

    Only query for sources above 1750 y pixels.

    ** And only query for sources about a certain MJD! **

    This is so messy, might only be able to make these on CR-corrected
    data.  Refer to `plot_cteslope_vs_flashlvl` below.

    Parameters:

    Returns:

    Outputs:
    """

    # Set the figure.
    pylab.figure(figsize=(12.5,8.5))

    all_fluxessubbkgrds = []
    all_mnbkgrds = []
    for key, val in zip(image_dict.keys(), image_dict.values()):
        # Each value of the dictionary is a list of lists,
        # [[fluxes], [bkgrds], [ypixs]]
        fluxes = np.array(val[0])
        print(" ")
        print("now on image {}".format(key))
        print("Fluxes: {}".format(fluxes))
        totbkgrds = np.array(val[2])
        if local_bkgrd:
            mnbkgrds = val[1]
            print("Mean backgrounds: {}".format(mnbkgrds))
            fluxes_sub_totbkgrds = fluxes-totbkgrds
            pylab.scatter(mnbkgrds, fluxes_sub_totbkgrds, s=25, marker='.', color='green') 
            for i in range(len(mnbkgrds)):
                all_fluxessubbkgrds.append(fluxes_sub_totbkgrds[i])
                all_mnbkgrds.append(mnbkgrds[i])
        else:
            # Use global bkgrd of image,
            global_bkgrd = query_for_globalbkgrd(targname, imagename=key)
            print("Global bkgrd: {}".format(global_bkgrd))
            if global_bkgrd != []:
                bkgrds = np.ones(len(fluxes)) * global_bkgrd[0]
                pylab.scatter(bkgrds, fluxes, s=25, marker='.', color='green')

    # Take an average of fluxes for each background bin and plot.
    bkgrd_bins = np.arange(1, 30)
    all_mnbkgrds=np.array(all_mnbkgrds)
    all_fluxessubbkgrds = np.array(all_fluxessubbkgrds)
    for bkgrd_bin in bkgrd_bins:
        bkgrds_bin1_indices = np.where(all_mnbkgrds <= bkgrd_bin+1.)[0]
        bkgrds_in_bin_low = all_mnbkgrds[bkgrds_bin1_indices]
        bkgrds_bin2_indices = np.where(bkgrds_in_bin_low >= bkgrd_bin)[0]
        bkgrds_in_bin = bkgrds_in_bin_low[bkgrds_bin2_indices]

        fluxes_in_bin = (all_fluxessubbkgrds[bkgrds_bin1_indices])[bkgrds_bin2_indices]
        print(fluxes_in_bin)

        sigmaclipped_fluxes = sigmaclip(fluxes_in_bin, high=4, low=4)
        mean_fluxes_in_bin = np.mean(sigmaclipped_fluxes[0])
        stdev_fluxes_in_bin = np.std(sigmaclipped_fluxes[0]) / 10.
        print("Sigma clipped mean of bkgrd bin {}-{} is {}".format(bkgrd_bin, bkgrd_bin+1, 
            mean_fluxes_in_bin))
        pylab.scatter( ((bkgrd_bin+1) + bkgrd_bin)/2., mean_fluxes_in_bin, 
            s=50, color='black', marker='s' )
        pylab.errorbar(((bkgrd_bin+1) + bkgrd_bin)/2., 
            mean_fluxes_in_bin, 
            yerr=stdev_fluxes_in_bin, 
            color='gray', 
            alpha=0.75, marker=None, ls='None')

    # Complete the plot and save. 
    pylab.xlabel('Local Sky Background [e-]', fontsize=22, weight='bold') 
    pylab.ylabel('Flux [e-], aperture 3 pxl radius, y > 1750 pixels', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    # Add line at 100 ei flux.
    pylab.axhline(y=100,linestyle=':', color='grey')
    # Add line at 12 e- background.
    pylab.axvline(x=12, linestyle='--', color='red')
    title = "{} {} proposal{} exp'{}' ctecorr{} chip{}".format(
        param_dict['targname'], param_dict['filter'], param_dict['proposid'],
        param_dict['exp_length'], param_dict['ctecorr'], param_dict['chip'])
    pylab.title(title, fontsize=16)
    pylab.xlim([-1, 30])
    pylab.ylim([-50, 20000])
    #pylab.legend(scatterpoints=1, loc='upper left')
    #pylab.text(56500, .00055, 'Standard Error: stddev / sqrt(# sources)')

    # Create outname for plot and text file.
    outname = 'fluxVSbkgrd_{}_{}_{}_{}_ctecorr{}_ch{}'.format(
        param_dict['targname'], param_dict['filter'], param_dict['proposid'], 
        param_dict['exp_length'], param_dict['ctecorr'], param_dict['chip'])
    pylab.savefig(os.path.join(outloc, '{}.png'.format(outname)), 
        bbox_inches='tight')

    return '{}.png'.format(outname)


#-------------------------------------------------------------------------------# 

def plot_cteslope_vs_flashlvl(slopes, flashlvls, standerrs, param_dict, 
    outloc='', plot_errorbars=False, xlims=[-1,120], ylims=[-0.1,0.4]):
    """ Plots CTE slope vs flashlvl for target-filter-exposure length-
    ctecorr-aperture-fluxbin combination.

    This will be the alternative to `plot_flux_vs_bkgrd` above.
    Not ideal because the flashlvl is not an exact measure of actual 
    background. But we may then see an interesting trend.  

    Parameters:
        slopes : dict of lists
            Keys of epoch. Values of lists of CTE slopes.
        flashlvls : dict of lists
            Keys of epoch. Values of lists of CTE flashlvls.
        standerrs :dict of lists
            Keys of epoch. Values of lists of standard errors 
        param_dict : dict
            Dictionary containing image's header values and other
            useful information that will be used to fill in database
            and label plots and data files.
        outloc : string
            Path to the output location. Working directory by default.
        plot_errorbars : {True, False}
            Set to True to plot errorbars.  False by default.
        xlims : list of floats/ints
            The lower and upper bounds of the x-axis (flash levels).
        ylims : list of floats/ints
            The lower and upper bounds of the y-axis (CTE slope). 

    Returns:
        String. Name of the plot. 

    Outputs:
        PNG plot.
        "cteVSflashlvl_<targname>_<filter>_<exposure length>_ctecorr<True/False>_r<aperture>_fluxbin<fluxbin>.png"

    """
   # Create outname for plot and text file.
    outname = 'cteVSflashlvl_{}_{}_{}_ctecorr{}_r{}_fluxbin{}'.format(
        param_dict['targname'], param_dict['filter'], 
        param_dict['exp_length'], param_dict['ctecorr'], 
        param_dict['aperture'], 
        param_dict['fluxbin'])

    # One color for each epoch, cycling using itertools.
    colors = itertools.cycle(['brown', 'red', 'darkgoldenrod',
                              'gold', 'lime', 'green', 
                              'cyan', 'blue', 'purple', 'magenta', 
                              'grey', 'black'])

    # Set plot.
    fig=pylab.figure(figsize=(12.5,8.5))

    # Loop over available epochs.
    for epoch in slopes.keys():

        if slopes[epoch] != []:
            # Change all Nones to NaNs.
            slopes_epoch = np.array(slopes[epoch], dtype=np.float)
            standerrs_epoch = np.array(standerrs[epoch], dtype=np.float)

            # Set the next color in sequence.
            color = next(colors)

            # Get x and y values into desired formats.
            flashlvls_epoch = np.array(flashlvls[epoch])
            # Correct the CTE measurements for the chip-height (2048 pixels).
            corrected_cteslopes = (slopes_epoch / 2.) * 2048.
            corrected_standerrs = (standerrs_epoch / 2.) * 2048.

            # Plot!
            pylab.scatter(flashlvls_epoch, corrected_cteslopes, 
                s=55, marker='D', alpha=0.65, color = color,
                label='MJD {}'.format(int(epoch)))
            if plot_errorbars:
                pylab.errorbar(flashlvls_epoch, corrected_cteslopes,
                    yerr=corrected_standerrs, 
                    color='gray', alpha=0.75, marker=None, ls='None')

    # Complete the plot and save. 
    pylab.xlabel('Header FLASHLVL Value [-e]', fontsize=22, weight='bold') 
    pylab.ylabel('CTE loss [flux / 2048 pxl]', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    # Add lines at major y tickmarks. Must define by hand...
    for ymaj in np.arange(11)*0.1:
        pylab.axhline(y=ymaj,linestyle='--',color='grey')
    # Add vertical line at 12 e- background.
    pylab.axvline(x=12., linestyle='--', color='red')
    title = "{} {} explen'{}' ctecorr{} ap{} fluxbin{}".format(
        param_dict['targname'], param_dict['filter'], 
        param_dict['exp_length'], 
        param_dict['ctecorr'], param_dict['aperture'],
        param_dict['fluxbin'])
    pylab.title(title, fontsize=16)
    if xlims != []:
        pylab.xlim(xlims)
    if ylims != []:
        pylab.ylim(ylims)

    try:
        pylab.legend(scatterpoints=1, loc='best', ncol=3)
    except:
        print("Unable to make legend for {}.png!".format(outname))

    pylab.savefig(os.path.join(outloc, '{}.png'.format(outname)), 
        bbox_inches='tight')

    return '{}.png'.format(outname)


#-------------------------------------------------------------------------------# 

def plot_180_slope_vs_expt(slopes_chip1, standerrs_chip1, expts_chip1, 
    slopes_chip2, standerrs_chip2, expts_chip2, param_dict, outloc='',
    chip1_invert=True, plot_errorbars=False):
    """Plots the CTE slopes vs exposure times of the full 180-degree dataset.
    The chips are on the same plot, color-coded. 

    Parameters:
        slopes_chip1 : dictionary
            Keys of filter name. Values of list of Chip 1 CTE slopes.
        standerrs_chip1 : dictionary
            Keys of filter name. Values of list of Chip 1 CTE slope 
            standard errors.
        expts_chip1 : dictionary
            Keys of filter name. Values of list of Chip 1 exposure times.
        slopes_chip1=2 : dictionary
            Keys of filter name. Values of list of Chip 2 CTE slopes.
        standerrs_chip2 : dictionary
            Keys of filter name. Values of list of Chip 2 CTE slope 
            standard errors.
        expts_chip2 : dictionary
            Keys of filter name. Values of list of Chip 2 exposure times.
        param_dict : dictionary
            Paramters of the dataset.
        outloc : string
            Path to the output location. Working directory by default.
        chip1_invert : {True, False}
            Multipy -1.0 to values of the chip1 fluxes. True by default.
            (Chip1 slopes negative because direction of readout in opposite 
            direction from chip2's direction.)
        plot_errorbars : {True, False}
            Set to True to plot errorbars.  False by default.

    Returns:
        Name of the plots.

    Outputs:
        PNG plot.
        "180cteVSexpt_ctecorr<True/False>_r<aperture>_<fluxbin>.png"
    """
    # Create outname for plot and text file.
    outname = '180cteVSexpt_ctecorr{}_r{}_{}'.format(
        param_dict['ctecorr'], param_dict['aperture'], param_dict['fluxbin'])

    # Set plot.
    fig=pylab.figure(figsize=(12.5,8.5))
 
    markers = ['o', '^']

    # Loop over filters and their markers.
    for filt, marker in zip(slopes_chip1.keys(), markers):
        print(filt)
        if slopes_chip1[filt] != [] and slopes_chip1[filt] != [[None]]:
            # Change all Nones to NaNs.
            if chip1_invert:
                slopes_chip1_filt = (-1.0)*(np.array(slopes_chip1[filt], dtype=np.float))
            else:
                slopes_chip1_filt = np.array(slopes_chip1[filt], dtype=np.float)         
            standerrs_chip1_filt = np.array(standerrs_chip1[filt], dtype=np.float)

            # Correct the CTE measurements for the chip-height (2048 pixels).
            corrected_cteslopes_chip1 = (slopes_chip1_filt / 2.) * 2048.
            corrected_standerrs_chip1 = (standerrs_chip1_filt / 2.) * 2048.

            # Plot!
            pylab.scatter(expts_chip1[filt], corrected_cteslopes_chip1, 
                s=90, marker=marker, alpha=0.75, color='blue',
                label='{} chip1'.format(filt))
            if plot_errorbars:
                pylab.errorbar(expts_chip1[filt], corrected_cteslopes_chip1, 
                    yerr=corrected_standerrs_chip1, 
                    color='gray', alpha=0.75, marker=None, ls='None')

        if slopes_chip2[filt] != [] and slopes_chip2[filt] != [[None]]:
            # Change all Nones to NaNs.
            slopes_chip2_filt = np.array(slopes_chip2[filt], dtype=np.float)
            standerrs_chip2_filt = np.array(standerrs_chip2[filt], dtype=np.float)

            # Get values into desired formats.
            corrected_cteslopes_chip2 = (slopes_chip2_filt / 2.) * 2048.
            corrected_standerrs_chip2 = (standerrs_chip2_filt / 2.) * 2048.

            # Plot!
            pylab.scatter(expts_chip2[filt], corrected_cteslopes_chip2, 
                s=90, marker=marker, alpha=0.75, color='orange',
                label='{} chip2'.format(filt))
            if plot_errorbars:
                pylab.errorbar(expts_chip2[filt], corrected_cteslopes_chip2, 
                    yerr=corrected_standerrs_chip2, 
                    color='gray', alpha=0.75, marker=None, ls='None')

    # Complete the plot and save. 
    pylab.xlabel('exposure time [s]', fontsize=22, weight='bold') 
    pylab.ylabel('CTE loss [flux / 2048 pxl]', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)

    title = "ctecorr{} r{} fluxbin{}".format(param_dict['ctecorr'], 
        param_dict['aperture'], param_dict['fluxbin'])
    pylab.title(title, fontsize=16)

    pylab.axhline(y=0., linestyle='--', color='blue')
    # Add lines at major y tickmarks. Must define by hand...
    for ymaj in np.arange(6)*0.1: 
        pylab.axhline(y=ymaj,linestyle='--',color='grey')
    for ymaj in np.arange(6)*(-0.1):
        pylab.axhline(y=ymaj,linestyle='--',color='grey')
    #pylab.xlim(xlims)
    pylab.ylim([-0.1, 0.4])

    try:
        pylab.legend(scatterpoints=1, loc='upper center', ncol=2)
    except:
        print("couldn't create legend")

    pylab.savefig(os.path.join(outloc, '{}.png'.format(outname)), 
        bbox_inches='tight')

    return '{}.png'.format(outname)


#-------------------------------------------------------------------------------# 
# Wrapper plot functions that query the database and loop over images.
#-------------------------------------------------------------------------------# 

def plot_fluxratio_vs_ypos_setup(targname, proposid, dateobs, filt, exptime, 
                   chinject, flashlvl, ctecorr, postarg1, aperture=3,
                   chip=0, outloc='', write_ratios_to_file=False, 
                   insert_into_db=True):
    """ 
    Queries for lists to be plotted and passes them to the plotting function.
    Creates for each image pair and each fluxbin, the ratio of the fluxes 
    (chip1/chip2) over chip2's y-position.

    Parameters:
        targname : string
            Name of the target cluster.
        proposid : string
            Proposal number. 
        dateobs : float
            MJD date of observation. 
        filt : string
            Name of the filter.
        exptime : float
            Expsosure time.
        chinject : string
            Charge injection.
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
        postarg1 : float
            Postarg1.
        aperture : int
            The aperture of photometry. (3 is nominal)
        chip : int
            If querying for 180-degree datasets, set this parameter for 
            the chip number of desired pair.
        outloc : string
            Path to desired output location.
        write_ratios_to_file : {True, False}
            Set to True to record the ratios for each flux bin in a file 
            '*_slopes.txt'.
        insert_into_db : {True, False}
            Set on to insert slopes into Results table for each fluxbin.

    Returns:
        nothing

    Outputs:
        For each image pair, for each fluxbin,

         PNG plot. If nominal dataset,
        "<imagename1>_<imagename2>_r<aperture>_slopes.png"
        And if 'write_ratios_to_file' set to True, text file,
        "<imagename1>_<imagename2>_r<aperture>_slopes.txt"

        If 180-degree dataset,
        "<imagename1>_<imagename2>_r<aperture>_ch<1/2>_slopes.png"    
        And if 'write_ratios_to_file' set to True, text file,
        "<imagename1>_<imagename2>_r<aperture>_ch<1/2>_slopes.txt"   

    """
    print("---")
    print("Querying FileInfo tables for targname={}, proposid={}, dateobs={}, \
          filter={}, exptime={}, chinject={}, flashlvl={}, ctecorr={}, \
          postarg1={}, aperture={}.".format(targname, proposid, dateobs, 
          filt, exptime, chinject, flashlvl, ctecorr, postarg1, aperture))
    print("---")

    # Query for all images that match the parameters. 
    # Use a different query if doing 180-degree dataset (NGC6583).
    if '6583' in targname:
        # Note that in this dataset, the two images are actually the SAME
        # chip. But to make the pairs work through rest of function, 
        # just pretend.
        imagename_chip1_ls, imagename_chip2_ls = query_for_180pair(
            targname=targname, filt=filt, exptime=exptime, 
            ctecorr=ctecorr, chip=chip)

    else:
        imagename_chip1_ls, imagename_chip2_ls = query_for_pair(
            targname=targname, proposid=proposid, dateobs=dateobs, 
            filt=filt, exptime=exptime, chinject=chinject, flashlvl=flashlvl, 
            ctecorr=ctecorr, postarg1=postarg1)
        print("Query returned the chip1 imagenames {} and chip2 imagenames {}."\
            .format(imagename_chip1_ls, imagename_chip2_ls))

    # Check that a pair was returned. 
    if (len(imagename_chip1_ls) == len(imagename_chip2_ls)) and \
       (imagename_chip1_ls != [] and imagename_chip2_ls != []):

        # Loop through image pairs. 
        # There should be only one image per chip for the nominal datasets.
        # There may be more than one for 180-degree dataset
        for imagename_chip1, imagename_chip2 in zip(imagename_chip1_ls, imagename_chip2_ls):

            fluxes_chip1, bkgrds_chip1, ypos_chip1, master_id_chip1 = \
                query_for_flux_by_imagename(targname, 
                                            imagename=imagename_chip1, 
                                            aperture=aperture, 
                                            bkgrd_returned='tot')
            fluxes_chip2, bkgrds_chip2, ypos_chip2, master_id_chip2  = \
                query_for_flux_by_imagename(targname, 
                                            imagename=imagename_chip2, 
                                            aperture=aperture, 
                                            bkgrd_returned='tot')

            # Check that the queries returned results and the same number for each chip.
            if (len(fluxes_chip1) == len(fluxes_chip1)) and (fluxes_chip1 != [] and fluxes_chip2 != []):

                # Make sure that the master_ids from both chips match. 
                matched_indices_chip1, matched_indices_chip2 = match_master_ids(master_id_chip1, master_id_chip2)

                master_id_chip1 = np.array(master_id_chip1)[matched_indices_chip1]
                master_id_chip2 = np.array(master_id_chip2)[matched_indices_chip2]

                # Pull out entries that have a matched master ID for both chips
                fluxes_chip1 = np.array(fluxes_chip1)[matched_indices_chip1]
                fluxes_chip2 = np.array(fluxes_chip2)[matched_indices_chip2]
                bkgrds_chip1 = np.array(bkgrds_chip1)[matched_indices_chip1]
                bkgrds_chip2 = np.array(bkgrds_chip2)[matched_indices_chip2]    
                ypos_chip1 = np.array(ypos_chip1)[matched_indices_chip1] # just in case want it later
                ypos_chip2 = np.array(ypos_chip2)[matched_indices_chip2]

                # Also input a param_dict so that can have a more informative title.
                param_dict = {'targname':targname, 'proposid':proposid, 
                              'dateobs':dateobs, 'filter':filt, 
                              'exptime':exptime, 'chinject':chinject, 
                              'flashlvl':flashlvl, 'ctecorr':ctecorr, 
                              'aperture':aperture}
                imagename_1, imagename_2, slope_file = plot_fluxratio_vs_ypos(
                    imagename_chip1=imagename_chip1, 
                    imagename_chip2=imagename_chip2, 
                    fluxes_chip1=fluxes_chip1, 
                    fluxes_chip2=fluxes_chip2, 
                    bkgrds_chip1=bkgrds_chip1, 
                    bkgrds_chip2=bkgrds_chip2, 
                    ypos_chip2=ypos_chip2, 
                    master_id_chip2=master_id_chip2, 
                    chip=chip,
                    param_dict=param_dict, 
                    write_ratios_to_file=write_ratios_to_file)

                # Update Results table in database.
                if insert_into_db:
                    results_dict = return_results_dict(imagename_1=imagename_1, 
                                                       imagename_2=imagename_2, 
                                                       aperture=aperture,
                                                       slope_file=slope_file, 
                                                       outloc=outloc)
                    FileInfo, Master, Phot, Results = return_tables(targname)
                    update_results_table(Results, results_dict)

            
            else:
                print("Entries missing in Phot table for chip1 image {} and chip2 image {}."\
                    .format(imagename_chip1, imagename_chip2))

    elif (len(imagename_chip1_ls) != len(imagename_chip2_ls)) and \
         (imagename_chip1_ls != [] and imagename_chip2_ls != []):
        print("Too many matches for the pair! chip1 image {} and chip2 image {}."\
            .format(imagename_chip1_ls, imagename_chip2_ls))

    else:
        print("No pair found")


#-------------------------------------------------------------------------------# 

def plot_cteslope_vs_time_setup(exp_length, filt, flashlvl, ctecorr, 
    aperture=3, chinject='NO', postarg1=0, xlims=[], ylims=[], outloc='',
    write_to_file=True):
    """ Performs queries in Results tables to create the flux-ratio vs
    yposition slope vs time plots for each filter, exposure length, flashlvl, 
    ctecorr, and aperture combination.

    Parameters:
        exp_length : string
            Length of exposure. Either 'l' or 's'.
        filt : string
            Name of the filter.
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
        aperture : int
            The aperture of photometry. (3 is nominal)
        chinject : string
            Charge injection.
        postarg1 : float
            Postarg1.
        xlims : list of floats
            The lower and upper bounds of the x-axis (MJD time). 
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope). 
        outloc : string
            Path to desired output location.
        write_to_file : {True, False}
            Set on to write the slope and time lists to a text file.
            True by default.

    Returns:
        nothing

    Outputs:
        For each combination of filter, exposure length, flashlvl, ctecorr,
        and aperture a PNG plot,
        "cteVStime_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>.png"
    """
    print("Starting CTEslope_vs_time for exp_length {}, filter {}, flashlvl {}, ctecorr {}"\
        .format(exp_length, filt, flashlvl, ctecorr))
    slopes = create_fluxbin_dict()
    standerrs = create_fluxbin_dict()
    times = create_fluxbin_dict()

    fluxbins = ['500-2000', '2000-8000', '8000-32000']

    # Loop over all targets.
    for targname in slopes.keys():
        FileInfo, Master, Phot, Results = return_tables(targname)

        # Loop over all proposals.
        for proposid in program_list:
            # Find valid exptimes.
            exptimes_raw = query_for_exptimes(targname, proposid, filt)
            exptimes_raw = np.array(list(set(exptimes_raw)))
            if exp_length == 'l':
                exptimes = exptimes_raw[np.where(exptimes_raw > 60.)]
            elif exp_length == 's':
                exptimes = exptimes_raw[np.where(exptimes_raw <= 60.)]
            print(exptimes)
            # Loop over all exposure times.
            for exptime in exptimes:
                # Find valid dateobss.
                dateobss = query_for_dateobss(targname, proposid, filt, exptime)
                dateobss = list(set(dateobss))
                # Loop over all valid dateobs.
                for dateobs in dateobss:

                    # Query the FileInfo table for imagenames.
                    imagename_chip1, imagename_chip2 = query_for_pair(
                               targname, proposid, dateobs, filt, exptime=exptime, 
                               chinject=chinject, flashlvl=flashlvl, 
                               ctecorr=ctecorr, postarg1=postarg1)
                    print(imagename_chip1)
                    print(imagename_chip2)
                    if imagename_chip1 != [] and imagename_chip2 !=[]:

                        # Query the Results for the imagenames, for each fluxbin.
                        for fluxbin in fluxbins:
                            print("fluxbin: {}".format(fluxbin))

                            slope, slopestdev, numpoints = query_results_for_slopes(
                                targname=targname, 
                                imagename_chip1=imagename_chip1[0], 
                                imagename_chip2=imagename_chip2[0], 
                                fluxbin=fluxbin, 
                                aperture=aperture)
                

                            if slope != []:
                                
                                # Append slope for the flux bin.
                                ((slopes[targname])[fluxbin]).append(slope[0])
                                # Get the standard error.
                                if slopestdev[0] != None:
                                    if not isinstance(numpoints[0], int):
                                        numpoint = int.from_bytes(numpoints[0], byteorder='little')
                                    else:
                                        numpoint = numpoints[0]

                                    if numpoint != 0:
                                        standerr = slopestdev[0] / numpoint 
                                    else:
                                        standerr = np.nan

                                else:
                                    standerr = np.nan
                                ((standerrs[targname])[fluxbin]).append(standerr)

                                print("slopes dict: {}".format(slopes[targname]))
                                print("standerrs dict: {}".format(standerrs[targname]))

                                # If query successful, append dateobs to list.
                                ((times[targname])[fluxbin]).append(dateobs)
                                print(times)

                    else:
                        print("No entries in FileInfo for proposid {}, exp_length {}, \
                         dateobs {}, imagename_chip1 {}, imagename_chip2{}"\
                         .format(proposid, exp_length, dateobs, imagename_chip1, \
                            imagename_chip2))
                        

    param_dict = {'filter':filt, 'exp_length':exp_length, 'flashlvl':flashlvl, 
                   'ctecorr':ctecorr, 'aperture':aperture}

    # Only if there are, in fact, results, create a plot.
    if not empty_tree(slopes['ngc104']) or not empty_tree(slopes['ngc6791']):
        plot_cteslope_vs_time(slopes=slopes, 
                              times=times, 
                              standerrs=standerrs, 
                              param_dict=param_dict, 
                              outloc=outloc,
                              write_to_file=write_to_file)


#-------------------------------------------------------------------------------# 

def plot_cteslope_vs_logflux_setup(targnames, exp_length, filt, flashlvl=0, 
    ctecorr=False, aperture=3, chinject='NO', postarg1=0, chip=0,
    ylims=[-0.1,0.6], clipoutliers=False, skipbins=1,
    plot_errorbars=True, outloc='', chip1_invert=True):
    """ Performs queries in Results tables to create the flux-ratio vs 
    y-position slope vs log10 flux plots for each target, filter, exposure 
    length, flashlvl, ctecorr, and aperture combination.

    Parameters:
        targnames : list of strings
            List of the names of the target clusters.
        exp_length : string
            Length of exposure. Either 'l' or 's'.
        filt : string
            Name of the filter.
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
        aperture : int
            The aperture of photometry. (3 is nominal)
        chinject : string
            Charge injection.
        postarg1 : float
            Postarg1.
        chip : int
            If querying for 180-degree datasets, set this parameter for 
            the chip number of desired pair.
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope). 
        clipoutliers : {True, False}
            Turn on to clip 2*stddev outliers from polynomial fit.
            Not very affective. Recomended to play with 'skipbins' in 
            :func:``plot_cteslope_vs_logflux`` instead.
        plot_errorbars : {True, False}
            Set to True to plot errorbars.
        outloc : string
            Path to desired output location.
        chip1_invert : {True, False}
            Multipy -1.0 to values of the chip1 fluxes. True by default.
            (180-degree dataset's chip1 slopes negative because direction of 
            readout in opposite direction from chip2's direction.)

    Returns:
        nothing

    Outputs:
        For each combination of filter, exposure length, flashlvl, ctecorr,
        and aperture:

        PNG plot. If nominal dataset,
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>.png"

        If 180-degree dataset,
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>_ch<1/2>.png"

        Text files containing polynomial fit's 9 coefficients for each
        epoch.
        "cteVSlogflux_<targname>_<filter>_<explength>_pf<flashlvl>_ctecorr<True/False>_r<aperture>_mjd<epoch>_coeffs.txt"

    """

    # The fluxbins_hi and _lo are not the values that will be averaged 
    # to be plotted, but are just used to query for the appropriate flux ranges.

    # Each target gets its own plot.
    for targname in targnames:
        FileInfo, Master, Phot, Results = return_tables(targname)

        print("Starting CTEslope_vs_logflux for target {} exp_length {}, filter {}, flashlvl {}, ctecorr {}"\
            .format(targname, exp_length, filt, flashlvl, ctecorr))

        # Find valid dateobss, ie, epochs. These values are actually
        # averaged epochs; all epochs that occur within 30 days of
        # each other are averaged together. 
        epochs = find_epochs(targname)
        print("Found epochs: {}".format(epochs))

        # Create dictionaries to store results to be plotted.
        slopes = create_epoch_dict(epochs)
        standerrs = create_epoch_dict(epochs)
        flux_means = create_epoch_dict(epochs)
        all_fluxes = create_epoch_dict(epochs)

        # Loop over all epochs.
        for epoch in epochs:
            print("Looping over epoch {}".format(epoch))
            # Find valid exptimes.
            # The query searches 30 days before and after the averaged epoch.
            exptimes_raw = query_for_exptimes(targname, proposid='', 
                filt=filt, dateobs=epoch)
            exptimes_raw = np.array(list(set(exptimes_raw)))
            if exp_length == 'l':
                exptimes = exptimes_raw[np.where(exptimes_raw > 60.)]
            elif exp_length == 's':
                exptimes = exptimes_raw[np.where(exptimes_raw <= 60.)]
            print("Found exptimes: {}".format(exptimes))

            for exptime in exptimes:
                print("Looping over exptime {}".format(exptime))
                # Now query for fluxes (ap 3, chip 1) in each bin and the 
                # corresponding imagenames. Find the actual high and low 
                # value for each bin and take the average. 
                # It is these values that will be plotted on a log scale. 
                for flux_lo, flux_hi in zip(fluxbins_lo, fluxbins_hi):
                    print("Looping over flux_lo {} and flux_hi {}".format(flux_lo, flux_hi))
                    if '6583' in targname: #180-degree dataset
                        fluxes_per_image_all, imagenames_all = query_for_flux_range(
                            targname, dateobs=epoch, 
                            exptime=exptime, flashlvl=flashlvl, 
                            ctecorr=ctecorr, filt=filt, chinject=chinject, 
                            postarg1=postarg1, chip=chip, flux_lo=flux_lo,
                            flux_hi=flux_hi, aperture=aperture)
                        imagenames_sorted = sorted(imagenames_all)
                        imagenames_1 = [im for im in imagenames_sorted if '10' in im]
                        imagenames_2 = [im for im in imagenames_sorted if '11' in im]
                        print("Found imagenames_1 {}".format(imagenames_1))
                        print("Found imagenames_2 {}".format(imagenames_2))
                        
                        matched_indices = match_indices(imagenames_1, imagenames_all)
                        fluxes_per_image = np.array(fluxes_per_image_all)[matched_indices]

                    else:
                        fluxes_per_image, imagenames_1 = query_for_flux_range(
                            targname, dateobs=epoch, 
                            exptime=exptime, flashlvl=flashlvl, 
                            ctecorr=ctecorr, filt=filt, chinject=chinject, 
                            postarg1=postarg1, chip=1, flux_lo=flux_lo,
                            flux_hi=flux_hi, aperture=aperture)
                        print("Found imagenames_1 {}".format(imagenames_1))
                        print("Length imagenames_1: {} and length fluxes_per_image: {}".format(len(imagenames_1), len(fluxes_per_image)))
                        # Find the matching imagename for chip 2.
                        imagenames_2 = []
                        for imagename_1 in imagenames_1:
                            imagename_2 = query_for_matching_imagename(targname, imagename_1)
                            imagenames_2.append(imagename_2)
                        print("Found imagenames_2 {}".format(imagenames_2))

                    # Query for slopes of each imagename pair. Insert slopes
                    # into the epoch dictionaries. 
                    for imagename_1, imagename_2, fluxes in \
                        zip(imagenames_1, imagenames_2, fluxes_per_image):
                        if imagename_1 != [] and imagename_2 != []:
                            ## Should background be selected before query made???

                            print("Looping over imagename_1 {} and imagename_2 {}"\
                                .format(imagename_1, imagename_2))
                            fluxes_sigmaclipped = sigmaclip(fluxes, high=4, low=4)[0]
                            print("len fluxes: {}".format(len(fluxes)))
                            print("len sigmaclipped: {}".format(len(fluxes_sigmaclipped)))
                            flux_mean = np.mean(fluxes_sigmaclipped)

                            fluxbin = '{}-{}'.format(flux_lo, flux_hi)
                            slope, slopestdev, numpoints = query_results_for_slopes(
                                targname=targname, 
                                imagename_chip1=imagename_1, 
                                imagename_chip2=imagename_2, 
                                fluxbin=fluxbin, 
                                aperture=aperture)

                            if slope != []:
                                print("slope: {}".format(slope[0]))
                                # Append slope for the flux bin.
                                slopes[epoch].append(slope[0])
                                # Get the standard error.
                                if slopestdev[0] != None:
                                    if not isinstance(numpoints[0], int):
                                        numpoint = int.from_bytes(numpoints[0], byteorder='little')
                                    else:
                                        numpoint = numpoints[0]
                                    standerr = slopestdev[0] / np.sqrt(numpoint)
                                else:
                                    standerr = np.nan
                                standerrs[epoch].append(standerr)

                                # If query successful, append the flux_mean
                                print("Flux mean: {}".format(flux_mean))
                                flux_means[epoch].append(flux_mean)
                                all_fluxes[epoch].append(fluxes_sigmaclipped)

                            else:
                                print("No entries in FileInfo for target {} exp_length {}, filter {}, flashlvl {}, ctecorr {}"
                                 .format(targname, exp_length, filt, flashlvl, ctecorr))
                            

        if '6791' in targname:
            # Low flux bins result in bad fit for this cluster.
            skipbins=2
        else:
            # The other two clusters get best fit with just first bin skipped.
            skipbins=1

        param_dict = {'targname':targname, 'filter':filt, 'exp_length':exp_length, 
                       'flashlvl':flashlvl, 'ctecorr':ctecorr, 'aperture':aperture}

        if not empty_tree(slopes):
            print(" ")
            print(slopes)
            print(flux_means)
            plot_cteslope_vs_logflux(
                targname=targname, 
                slopes=slopes, 
                flux_means=flux_means, 
                all_fluxes=all_fluxes, 
                standerrs=standerrs, 
                param_dict=param_dict, 
                chip=chip,
                outloc=outloc, 
                skipbins=skipbins,
                clipoutliers=clipoutliers,
                plot_errorbars=plot_errorbars,
                chip1_invert=chip1_invert)


#-------------------------------------------------------------------------------# 

def plot_flux_vs_bkgrd_setup(exp_length, ctecorr=False, local_bkgrd=True,
    xlims=[], ylims=[], outloc=''):
    """ This should emulate Figure 3 of Jay's ISR. 

    Only query for sources above 1750 y pixels.

    Do a pre-query on long exposures for stars with flux ~= 100.
    Take the master ids and query on them on the short exposures?
    Do this only on ngc104.

    Actually doesn't work too well (should resurrect after correct 
        images for cosmic rays...)
    See `plot_cteslope_vs_flashlvl_setup`

    Parameters:
        exp_length : string
            Length of exposure. Either 'l' or 's'.
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        local_bkgrd : {True, False}

        xlims : list of floats
            The lower and upper bounds of the x-axis (bkgrd). 
        ylims : list of floats
            The lower and upper bounds of the y-axis (flux). 
        outloc : string
            Path to desired output location.
    
    Returns:
        nothing

    Outputs:


    """
    
    # Query for all sources with ypos > 1750 pixels.
    # The query goes over all postflash levels.
    for targname in target_list[:1]: 
        for proposid in program_list:
            for filt in filter_list:
                for chip in [1,2]:
                    exptimes_raw = query_for_exptimes(targname=targname, 
                        proposid=proposid, filt=filt)
                    exptimes_raw = np.array(list(set(exptimes_raw)))
                    if exp_length == 'l':
                        exptimes = exptimes_raw[np.where(exptimes_raw > 60.)]
                    elif exp_length == 's':
                        exptimes = exptimes_raw[np.where(exptimes_raw <= 60.)]

                    for exptime in exptimes:
                        image_dict = query_for_fluxes_bkgrds_by_ypos(
                            targname, proposid=proposid, filt=filt, 
                            ctecorr=ctecorr, 
                            exptime=exptime, chip=chip,
                            chinject='NO', postarg1=0, ypos=1750)

                        param_dict = {'targname' : targname, 
                            'proposid' : proposid, 'chip' : chip, 
                            'filter' : filt, 'ctecorr' : ctecorr, 
                            'exp_length' : exp_length, 'exptime' : exptime}
                        plot_flux_vs_bkgrd(image_dict, param_dict, local_bkgrd)


#-------------------------------------------------------------------------------# 

def plot_cteslope_vs_flashlvl_setup(exp_length, filt, ctecorr, 
    aperture=3, chinject='NO', postarg1=0, xlims=[-1,120], 
    ylims=[-0.1,0.3], outloc=''):
    """ Creates plot of CTE-slope vs flashlvl for each target, filter,
    exposure length, ctecorr, and aperture combination.

    Parameters:
        exp_length : string
            Length of exposure. Either 'l' or 's'.
        filt : string
            Name of the filter.
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        aperture : int
            Size of the aperture in pixels.
        chinject : string
            Charge injection.
        postarg1 : float
            Postarg1.
        xlims : list of floats/ints
            The lower and upper bounds of the x-axis (flashlvl). 
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope). 
        outloc : string
            Path to the output location. Working directory by default.

    Returns:
        nothing

    Outputs:
       PNG plot.
        "cteVSflashlvl_<targname>_<filter>_<exposure length>_ctecorr<True/False>_r<aperture>_fluxbin<fluxbin>.png"
    """

    # Each target gets its own plot.
    for targname in target_list[:2]:
        FileInfo, Master, Phot, Results = return_tables(targname)

        print("Starting CTEslope_vs_logflux for target {} exp_length {}, filter {}, ctecorr {}"\
            .format(targname, exp_length, filt, ctecorr))

        # Find valid dateobss, ie, epochs. These values are actually
        # averaged epochs; all epochs that occur within 30 days of
        # each other are averaged together. 
        epochs = find_epochs(targname)
        print("Found epochs: {}".format(epochs))

        # Loop over flux bins.
        for flux_lo, flux_hi in zip(fluxbins_lo, fluxbins_hi):
            fluxbin = '{}-{}'.format(flux_lo, flux_hi)
            print("Looping over flux_lo {} and flux_hi {}"\
                .format(flux_lo, flux_hi))

            # Create dictionaries to store results to be plotted.
            slopes = create_epoch_dict(epochs)
            standerrs = create_epoch_dict(epochs)
            flashlvls = create_epoch_dict(epochs)

            # Loop over all epochs.
            for epoch in epochs:
                print("Looping over epoch {}".format(epoch))
                # Find valid exptimes.
                # The query searches 30 days before and after the averaged epoch.
                exptimes_raw = query_for_exptimes(targname, proposid='', 
                    filt=filt, dateobs=epoch)
                exptimes_raw = np.array(list(set(exptimes_raw)))
                if exp_length == 'l':
                    exptimes = exptimes_raw[np.where(exptimes_raw > 60.)]
                elif exp_length == 's':
                    exptimes = exptimes_raw[np.where(exptimes_raw <= 60.)]
                print("Found exptimes: {}".format(exptimes))

                # Loop over exposure times.
                for exptime in exptimes:
                    print("Looping over exptime {}".format(exptime))

                    # Loop over flash levels.
                    for flashlvl in flashlvl_list:
                        # remove the flux_lo and flux_hi sorting???
                        print("Looping over flashlvl: {}".format(flashlvl))

                        imagenames_1, imagenames_2 = query_for_pair(
                            targname, proposid='', dateobs=epoch, filt=filt, 
                            exptime=exptime, chinject=chinject, flashlvl=flashlvl, 
                            ctecorr=ctecorr, postarg1=postarg1)
                        print("Found imagenames_1 {}".format(imagenames_1))
                        print("Found imagenames_2 {}".format(imagenames_2))

                        # Query for slopes of each imagename pair. Insert slopes
                        # into the epoch dictionaries. 
                        for imagename_1, imagename_2 in zip(imagenames_1, imagenames_2):
                            if imagename_1 != [] and imagename_2 != []:
                                print("Looping over imagename_1 {} and imagename_2 {}"\
                                    .format(imagename_1, imagename_2))

                                slope, slopestdev, numpoints = query_results_for_slopes(
                                    targname=targname, 
                                    imagename_chip1=imagename_1, 
                                    imagename_chip2=imagename_2, 
                                    fluxbin=fluxbin, 
                                    aperture=aperture)

                                if slope != []:
                                    print("slope: {}".format(slope[0]))
                                    # Append slope for the flux bin.
                                    slopes[epoch].append(slope[0])
                                    # Get the standard error.
                                    if slopestdev[0] != None:
                                        standerr = slopestdev[0] / np.sqrt(numpoints[0])
                                    else:
                                        standerr = np.nan
                                    standerrs[epoch].append(standerr)

                                    # If query successful, append the flashlvl
                                    flashlvls[epoch].append(float(flashlvl))

                                else:
                                    print("No entries in FileInfo for target {} exp_length {}, filter {}, flashlvl {}, ctecorr {}"
                                     .format(targname, exp_length, filt, 
                                             flashlvl, ctecorr))
                            else:
                                print("Unable to loop over imagename_1 {} and imagename_2 {}"\
                                    .format(imagename_1, imagename_2))
                      

            # Now create plot for each targname and fluxbin 
            param_dict = {'targname':targname, 'filter':filt, 'exp_length':exp_length, 
                           'ctecorr':ctecorr, 'aperture':aperture, 'fluxbin':fluxbin}
                           #need to remove fluxbin sorting?

            # Create a plot only if query returned data.
            if not empty_tree(slopes):
                print(" ")
                print("slopes: {}".format(slopes))
                print("flashlvls: {}".format(flashlvls))
                if targname == 'ngc6791' or exp_length == 's':
                    xlims=[-1,13]
                plot_cteslope_vs_flashlvl(slopes=slopes, flashlvls=flashlvls, 
                    standerrs=standerrs, param_dict=param_dict, outloc=outloc,
                    plot_errorbars=True,
                    xlims=xlims, ylims=ylims)


#-------------------------------------------------------------------------------# 

def plot_180_slope_vs_expt_setup(ctecorr, aprads=[3,5], outloc='', 
    chip1_invert=True):
    """

    Parameters:
        ctecorr : {True, False}
            True for FLCs (pixel-based CTE-corrected)
            False for FLTs (not pixel-based CTE-corrected)
        aprads : list of ints/floats
            The aperture radii (pixels) on which to perform photometry.
        outloc : string
            Path to the output location. Working directory by default.
        chip1_invert : {True, False}
            Multipy -1.0 to values of the chip1 fluxes. True by default.
            (Chip1 slopes negative because direction of readout in opposite 
            direction from chip2's direction.)

    Returns:
        nothing.

    Outputs:
        For each combination of target, filter, ctecorr, aperture, and fluxbin, 
        a PNG plot.

        "180cteVSexpt_ctecorr<True/False>_r<aperture>_<fluxbin>.png"    
    """
    targname = 'ngc6583'

    for aperture in aprads:
        for fluxbin in ['500-2000', '2000-8000', '8000-32000']:
            print("fluxbin: {}".format(fluxbin))
            param_dict = {'aperture' : aperture, 'ctecorr' : ctecorr, 
                'fluxbin' : fluxbin}

            slopes_chip1 = make_OrderedDict({'F606W':[], 'F502N':[]})
            standerrs_chip1 = make_OrderedDict({'F606W':[], 'F502N':[]})
            expts_chip1 = make_OrderedDict({'F606W':[], 'F502N':[]})
            slopes_chip2 = make_OrderedDict({'F606W':[], 'F502N':[]})
            standerrs_chip2 = make_OrderedDict({'F606W':[], 'F502N':[]})
            expts_chip2 = make_OrderedDict({'F606W':[], 'F502N':[]})

            for filt in ['F606W', 'F502N']:
                print("filter: {}".format(filt))
                for exptime in [60.0, 348.0]:
                    print("exptime: {}".format(exptime))

                    imagenames_1_chip1, imagenames_2_chip1 = query_for_180pair(
                        targname=targname, 
                        filt=filt, 
                        exptime=exptime, 
                        ctecorr=ctecorr, 
                        chip=1)
                    
                    if imagenames_1_chip1 != [] and imagenames_2_chip1 != []:
                        for imagename_1_chip1, imagename_2_chip1 in \
                        zip(imagenames_1_chip1, imagenames_2_chip1):
                            slope_chip1, slopestdev_chip1, numpoints_chip1 = \
                            query_results_for_slopes(
                                targname=targname, 
                                imagename_chip1=imagename_1_chip1, 
                                imagename_chip2=imagename_2_chip1, 
                                fluxbin=fluxbin, 
                                aperture=aperture)
                            
                            slopes_chip1[filt].append(slope_chip1)
                            standerrs_chip1[filt].append(slopestdev_chip1)
                            expts_chip1[filt].append(exptime)    


                    imagenames_1_chip2, imagenames_2_chip2 = query_for_180pair(
                        targname=targname, 
                        filt=filt, 
                        exptime=exptime, 
                        ctecorr=ctecorr, 
                        chip=2)

                    if imagenames_1_chip1 != [] and imagenames_2_chip1 != []:
                        for imagename_1_chip2, imagename_2_chip2 in \
                        zip(imagenames_1_chip2, imagenames_2_chip2):
                            slope_chip2, slopestdev_chip2, numpoints_chip2 = \
                            query_results_for_slopes(
                                targname=targname, 
                                imagename_chip1=imagename_1_chip2, 
                                imagename_chip2=imagename_2_chip2, 
                                fluxbin=fluxbin, 
                                aperture=aperture)
                            
                            slopes_chip2[filt].append(slope_chip2)
                            standerrs_chip2[filt].append(slopestdev_chip2)
                            expts_chip2[filt].append(exptime)                                         

            print(slopes_chip1)
            print(slopes_chip2)
            if (not empty_tree(slopes_chip1)) and (not empty_tree(slopes_chip2)):
                plot_180_slope_vs_expt(
                    slopes_chip1=slopes_chip1, 
                    standerrs_chip1=standerrs_chip1, 
                    expts_chip1=expts_chip1, 
                    slopes_chip2=slopes_chip2, 
                    standerrs_chip2=standerrs_chip2, 
                    expts_chip2=expts_chip2, 
                    param_dict=param_dict, 
                    outloc=outloc)


#-------------------------------------------------------------------------------# 
# Helper functions.
#-------------------------------------------------------------------------------# 

def make_OrderedDict(dictionary):
    """ Makes an OrderedDict from given dictionary already containing
    keys and values.

    Parameters:
        dictionary : dictionary
            Any kind.

    Returns:
        ordered_dict : OrderedDict
            The ordered dict version of the input dictionary.

    Outputs:
        nothing.
    """
    ordered_dict = OrderedDict()
    for key, val in zip(dictionary.keys(), dictionary.values()):
        ordered_dict[key] = val

    return ordered_dict


#-------------------------------------------------------------------------------# 

def match_master_ids(master_id_chip1, master_id_chip2):
    """
    Returns indices for both chips which contain only matches.
    The returned indices should be the same size. 

    Parameters:
        master_id_chip1 : list or array
            Master IDs for chip1.
        master_id_chip2 : list or array
            List of master IDs for chip2.

    Returns:
        matched_indices_chip1 : list
            Indices from 'master_id_chip1', whose entries have matches in 
            'master_id_chip2'.
        matched_indices_chip2 : list
            Indices from 'master_id_chip2', whose entries have matches in 
            'master_id_chip1'.

    Outputs:
        nothing.
    """

    # Find the intersection of the two lists of master IDs.
    master_id_intersection = list(set(master_id_chip1) & set(master_id_chip2))

    # Get the indices of the matched master IDs from each chip.
    matched_indices_chip1 = match_indices(master_id_intersection, master_id_chip1) 
         #[i for i,x in enumerate(master_id_chip1) if x in master_id_intersection]
    matched_indices_chip2 = match_indices(master_id_intersection, master_id_chip2) 
        #[i for i,x in enumerate(master_id_chip2) if x in master_id_intersection]

    return matched_indices_chip1, matched_indices_chip2


#-------------------------------------------------------------------------------# 

def match_indices(shortened_list, primary_list):
    """
    Returns the 'primary_list's indices that correspond to the matching
    values in the 'shortened_list'. Assumes all values are unique.
    (For use in the External CTE monitor, matches are between RA and Decs
    and for a given source we can assume uniqueness.)

    Parameters:
        shortened_list : list of anything
            A shortened 'primary_list' of whose values you want the 
            'primary_list's indices. Assumes all values are unique.
        primary_list : list of anything 
            The original list of which the 'shortened_list' is a subset.
            Assumes all values are unique.
    Returns:
        matched_indices : list of ints 
            The 'primary_list's indices that correspond to the matching
            values in the 'shortened_list'.

    Outputs:
        nothing
    """
    matched_indices = [i for i,x in enumerate(primary_list) if x in shortened_list]
    return matched_indices


#-------------------------------------------------------------------------------# 

def create_fluxbin_dict():
    """ Creates an ordered dictionary of ordered dictionaries.

    Parameters:
        nothing
        
    Returns:
        fluxbin_dict : dict of dicts
            First keys of target. Second keys of flux bin. Values of 
            empty lists.

    Outputs:
        nothing
    """
    fluxbin_dict = OrderedDict()

    fluxbin_dict['ngc104'] = OrderedDict()
    (fluxbin_dict['ngc104'])['500-2000'] = []
    (fluxbin_dict['ngc104'])['2000-8000'] = []
    (fluxbin_dict['ngc104'])['8000-32000'] = []

    fluxbin_dict['ngc6791'] = OrderedDict()
    (fluxbin_dict['ngc6791'])['500-2000'] = []
    (fluxbin_dict['ngc6791'])['2000-8000'] = []
    (fluxbin_dict['ngc6791'])['8000-32000'] = []

    return fluxbin_dict

#-------------------------------------------------------------------------------# 

def create_epoch_dict(epochs):
    """ Creates an ordered dictionary of epoch keys and a list of lists,
    where the first list is the average flux of the bin, and second list
    the CTE loss slope.

    Parameters:
        epochs : list of ints/floats        
        
    Returns:
        epoch_dict : dict of lists
            Keys of epoch. Values of empty lists.

    Outputs:
        nothing
    """
    epoch_dict = OrderedDict()

    for epoch in epochs:

        epoch_dict[epoch] = []

    return epoch_dict

#-------------------------------------------------------------------------------# 

def empty_tree(input_dict):
    """Recursively iterate through values in nested dictionaries of lists.
    Returns True if the bottom level dictionaries' items are all
    empy lists.  Returns False if items are found.

    References:
http://stackoverflow.com/questions/1593564/python-how-to-check-if-a-nested-list-is-essentially-empty

    Parameters:
        input_dict : dict of lists
            Can be dict of dicts of lists. And so on.

    Returns:
        {True, False}

    Outputs:
        nothing
    """
    for item in input_dict.values():
        if isinstance(item, dict):
            empty_tree(input_dict)
        elif item != []:
            return False
    return True


#-------------------------------------------------------------------------------# 

def find_epochs(targname):
    """
    Use 'epochs' to query for fluxes?
    Then 'combined_epochs' to make sure all query results get placed 
    in right epoch lists.

    Return list of combined epochs, which can use to query
    each at a time for imagenames.

    Parameters:
        targname : string
            Name of the target cluster.

    Returns:
        combined_epochs : 

    Outputs:
        nothing
    """
    epochs = query_for_all_dateobss(targname)
    epochs = sorted(list(set(epochs)))
    print(epochs)
    # Average together those that fall within same month.
    combined_epochs = []
    i = 0
    while i <= len(epochs)-1:
        if i == len(epochs)-1:
            # Check the last value.
            combined_epochs.append(epochs[i])
            i += 1
        elif ((epochs[i+1] - epochs[i]) < 30):
            average_epoch = round(np.mean([epochs[i+1], epochs[i]]))
            combined_epochs.append(average_epoch)
            i += 2
        else:
            combined_epochs.append(epochs[i])
            i += 1
    print(combined_epochs)

    return combined_epochs

#-------------------------------------------------------------------------------# 

def polyfit2d(x, y, z, order=2):
    """Performs 2-D polyfit by default (can be adjusted to other orders).
    Returns slope of the fit.

    Could become more sophisticated, especially in way it handles end-points.

    Parameters:
        x : list of floats/ints
            The x-values.
        y : list of floats/ints
            The y-values.
        z : 

        order : int
            Order of the polynomial fit. 2 by default.

    Returns:
        m : float
            The slope to the polynomial fit.

    Outputs:
        nothing

    References:
    http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent    
    """
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)

    return m

#-------------------------------------------------------------------------------# 

def polyval2d(x, y, m):
    """ Just something I was playing with. Feel free to ignore...

    Parameters:
        x :
        y :
        m :

    Returns:
        z :

    Outputs:
        nothing
    """
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


def polyfit2d_2(X, Y, Z):
    """ Another something I was playing with. Also feel free to ignore...
    http://stackoverflow.com/questions/33964913/equivalent-of-polyfit-for-a-2d-polynomial-in-python
    """
    print("**Using polyfit2d_2!!**")   

    X = X.flatten()
    Y = Y.flatten()

    A = np.array([X*0+1, X, Y, X**2, X**2*Y, X**2*Y**2, Y**2, X*Y**2, X*Y]).T
    B = Z.flatten()

    coeffs, r, rank, s = np.linalg.lstsq(A, B)
    
    return coeffs


#-------------------------------------------------------------------------------# 
# The Main. Use for testing.
#-------------------------------------------------------------------------------# 

if __name__=="__main__":

    targname='ngc6791'
    proposid='12348'
    dateobs=55466
    filt='F502N'
    exptime=360
    chinject='NO'
    flashlvl=0
    ctecorr=False
    postarg1=0
    aperture=3
    exp_length='s'

    #plot_cteslope_vs_flashlvl_setup(exp_length=exp_length, filt=filt, ctecorr=ctecorr, 
    #   aperture=3, chinject='NO', postarg1=0)

    plot_cteslope_vs_logflux_setup([targname], exp_length, filt, flashlvl, ctecorr, 
        aperture=3, chinject='NO', postarg1=0,
        clipoutliers=False,
        outloc='')

    #plot_cteslope_vs_time_setup(exp_length, filt, flashlvl, ctecorr, aperture=3,
    #    chinject='NO', postarg1=0)

    #plot_fluxratio_vs_ypos_setup(targname=targname, proposid=proposid, 
    #                       dateobs=dateobs, filt=filt, exptime=exptime, 
    #                       chinject=chinject, flashlvl=flashlvl, ctecorr=ctecorr, 
    #                       postarg1=postarg1, aperture=aperture, 
    #                       write_ratios_to_file=True, insert_into_db=False)
    #
    #plot_flux_vs_bkgrd_setup(exp_length='l', ctecorr=False, local_bkgrd=True)
