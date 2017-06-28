#! /usr/bin/env python

"""
This was the notebook 'cte180test_plots.ipynb'.  But it decided to die.
never rely on the notebook.

So this creates the plots of the logflux10 vs CTE degredation, overplotting
the 180-degree results onto the results from the nominal monitor from the 
same proposal. 

They can be found in the 2016 ISR by Gosmeyer, C. and Bagget, S.   
"""

# Modify this.
outloc = ''

from uvis_external_cte_plots import fit_empirical_model 
from uvis_external_cte_plots import create_epoch_dict
from uvis_external_cte_plots import empty_tree
from uvis_external_cte_plots import match_indices
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

def plot_cteslope_vs_logflux(targnames, slopes, flux_means, all_fluxes, 
    standerrs, param_dict, chips, outloc='', ylims=[-0.1,0.4],
    skipbins=1, plot_errorbars=False, chip1_invert=True):
    """ 
    Ratio flux slopes (CTE) over the log of flux (nominally at 
    aperture 3 radius).  Fits a 2nd degree polynomial to each epoch
    and returns the 9 coefficients to the fit.

    Parameters:
        targnames : list
            Name of the target cluster.
        slopes : dict of lists
            Keys of epoch. Values of lists of CTE slopes.
        flux_means : dict of lists
            Keys of epoch. Values of lists of mean fluxes.
        all_fluxes : dict of lists
            Keys of epoch. Values of lists of fluxes.
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
        skipfirstbin : {True, False}
            Turn on to skip the first fluxbin when performing polynomial fit.
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
    print(targnames)
    outname = 'cteVSlogflux_180test_{}_{}_pf{}_ctecorr{}_r{}'.format(
        param_dict['filter'], param_dict['exp_length'], 
        param_dict['flashlvl'], param_dict['ctecorr'],
        param_dict['aperture'])

    # one for each epoch
    colors_ls = ['#529DB7', '#7DB874', '#D92120', '#404096', '#E39C37']
    #colors_ls = ['magenta', 'orange', 'red', 'blue', 'green']
    colors = itertools.cycle(colors_ls)

    # Set plot.
    fig=pylab.figure(figsize=(12.5,8.5))
    
    # Set a count, after which start plotting lines in different 
    # style, so that no two epochs have same color and style.
    epoch_count = 0

    for epoch in slopes.keys():
        if targnames[epoch] != []:
            targname = targnames[epoch][0]
        else:
            targname = ''
        print('targname: {}'.format(targname))
        
        # For each epoch and flux bin, need to average all slopes.
        # Take the log10 of the average flux of the flux bin.
        if flux_means[epoch] != []:
            chip = chips[epoch][0]    
            if chip == 0:
                linestyle = '--'
            else:
                linestyle = '-'
            #if epoch_count < len(colors_ls):
            #    linestyle = '--'
            #else:
            #    linestyle = '-'
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

            # overplot the empirical model fit.
            model_fit = fit_empirical_model(log10fluxes=log10fluxes, 
                                            corrected_cteslopes=corrected_cteslopes, 
                                            epoch=epoch,
                                            outname=outname,
                                            outloc=outloc,
                                            skipbins=skipbins,
                                            write_fit_values=True,
                                            chip=chip)
            pylab.plot(model_fit[0], model_fit[1], color=color, 
                linewidth=2.5, linestyle=linestyle, label='MJD {}, {}, {}'.format(int(epoch), targname, chip))


    # Complete the plot and save. 
    pylab.xlabel('LOG Flux (ap 3 pxl radius) [e-]', fontsize=22,
        weight='bold') 
    pylab.ylabel('CTE loss [flux / 2048 pxl]', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    # Add lines at major y tickmarks. Must define by hand...
    for ymaj in np.arange(11)*0.1:
        pylab.axhline(y=ymaj,linestyle='--',color='grey')

    title = "{} explen'{}' pf{} ctecorr{} ap{}".format(
        param_dict['filter'], 
        param_dict['exp_length'], param_dict['flashlvl'], 
        param_dict['ctecorr'], param_dict['aperture'])        
    pylab.title(title, fontsize=16)
    pylab.xlim([2.5, 4.5])
    pylab.ylim(ylims)
    pylab.legend(scatterpoints=1, loc='upper right', ncol=2)

    pylab.savefig(os.path.join(outloc, '{}.png'.format(outname)), 
        bbox_inches='tight')

    return '{}.png'.format(outname)

#-------------------------------------------------------------------------------# 

def plot_cteslope_vs_logflux_setup(targnames, exp_length, filt, flashlvl=0, 
    ctecorr=False, aperture=3, chinject='NO', postarg1=0,
    ylims=[-0.1,0.4], skipbins=1, plot_errorbars=True, outloc='',
    chip1_invert=True):
    """ Performs queries in Results tables to create the flux-ratio vs 
    yposition slope vs time plot.

    The flux is from Chip 2. 
    There is NO ctecorrection.
    The flashlvl must be 0.

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
        skipfirstbin : {True, False}
            Turn on to skip the first fluxbin when performing polynomial fit.
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

    """

    # The fluxbins_hi and _lo are not the values that will be averaged 
    # to be plotted, but are just used to query for the appropriate flux ranges. 

    epochs = [55854.0, 55855.0, 56001.0, 56008.0, 56103.01, 56103.02, 56124.0, 56131.0]

    # Create dictionaries to store results to be plotted.
    slopes = create_epoch_dict(epochs)
    standerrs = create_epoch_dict(epochs)
    flux_means = create_epoch_dict(epochs)
    all_fluxes = create_epoch_dict(epochs)
    targnames_per_epoch = create_epoch_dict(epochs)
    chips = create_epoch_dict(epochs)    
    
    
    # Each target gets its own plot.
    for targname in targnames:
        FileInfo, Master, Phot, Results = return_tables(targname)

        print("Starting CTEslope_vs_logflux for target {} exp_length {}, filter {}, flashlvl {}, ctecorr {}"\
            .format(targname, exp_length, filt, flashlvl, ctecorr))

        # Find valid dateobss, ie, epochs. These values are actually
        # averaged epochs; all epochs that occur within 30 days of
        # each other are averaged together. 
            
        if targname == 'NGC104':
            epochs = [55855.0, 56008.0, 56131.0]
        elif targname == 'NGC6791':
            epoch = [55854.0, 56001.0, 56124.0]
        elif targname == 'NGC6583':
            epochs = [56103.01, 56103.02]        
        
        # Loop over all epochs.
        for epoch in epochs:
            print("Looping over epoch {}".format(epoch))
            # Find valid exptimes.
            # The query searches 30 days before and after the averaged epoch.
            exptimes_raw = query_for_exptimes(targname, proposid='12692', 
                filt=filt, dateobs=int(epoch))
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
                        if '.01' in str(epoch):
                            chip = 1
                        elif '.02' in str(epoch):
                            chip = 2
                        
                        fluxes_per_image_all, imagenames_all = query_for_flux_range(
                            targname, dateobs=int(epoch), 
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

                        print("Chip number: {}".format(chip))
                        # Query for slopes of each imagename pair. Insert slopes
                        # into the epoch dictionaries.
                        for imagename_1, imagename_2, fluxes in \
                            zip(imagenames_1, imagenames_2, fluxes_per_image):
                            if imagename_1 != [] and imagename_2 != []:
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
                                        standerr = slopestdev[0] / np.sqrt(numpoints[0])
                                    else:
                                        standerr = np.nan
                                    standerrs[epoch].append(standerr)

                                    # If query successful, append the flux_mean
                                    print("Flux mean: {}".format(flux_mean))
                                    flux_means[epoch].append(flux_mean)
                                    all_fluxes[epoch].append(fluxes)

                                    targnames_per_epoch[epoch].append(targname)
                                    chips[epoch].append(chip)

                                else:
                                    print("No entries in FileInfo for target {} exp_length {}, filter {}, flashlvl {}, ctecorr {}"
                                     .format(targname, exp_length, filt, flashlvl, ctecorr))
                                print(" - - - ")

                    else:
                        chip=0
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

                        print("Chip number: {}".format(chip))
                        # Query for slopes of each imagename pair. Insert slopes
                        # into the epoch dictionaries.
                        for imagename_1, imagename_2, fluxes in \
                            zip(imagenames_1, imagenames_2, fluxes_per_image):
                            if imagename_1 != [] and imagename_2 != []:
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
                                        standerr = slopestdev[0] / np.sqrt(numpoints[0])
                                    else:
                                        standerr = np.nan
                                    standerrs[epoch].append(standerr)

                                    # If query successful, append the flux_mean
                                    print("Flux mean: {}".format(flux_mean))
                                    flux_means[epoch].append(flux_mean)
                                    all_fluxes[epoch].append(fluxes)

                                    targnames_per_epoch[epoch].append(targname)
                                    chips[epoch].append(chip)

                                else:
                                    print("No entries in FileInfo for target {} exp_length {}, filter {}, flashlvl {}, ctecorr {}"
                                     .format(targname, exp_length, filt, flashlvl, ctecorr))
                                print(" - - - ")
                            

    param_dict = {'filter':filt, 'exp_length':exp_length, 
                   'flashlvl':flashlvl, 'ctecorr':ctecorr, 'aperture':aperture}

    if not empty_tree(slopes):
        print(" ")
        print(slopes)
        print(flux_means)
        print(targnames_per_epoch)
        print(chips)
        plot_cteslope_vs_logflux(
            targnames=targnames_per_epoch, 
            slopes=slopes, 
            flux_means=flux_means, 
            all_fluxes=all_fluxes, 
            standerrs=standerrs, 
            param_dict=param_dict, 
            chips=chips,
            outloc=outloc, 
            skipbins=skipbins, 
            plot_errorbars=plot_errorbars,
            chip1_invert=chip1_invert)


#-------------------------------------------------------------------------------# 
# The Main.
#-------------------------------------------------------------------------------# 

if __name__=='__main__':

    targnames = ['NGC104', 'NGC6583']

    for ctecorr in [True, False]:
        for exp_length in ['s', 'l']:
            for filt in ['F606W', 'F502N']:
                plot_cteslope_vs_logflux_setup(targnames, exp_length, filt, flashlvl=0, 
                    ctecorr=ctecorr, aperture=5, chinject='NO', postarg1=0,
                    ylims=[-0.1,0.4], skipbins=1, plot_errorbars=True, outloc=outloc,
                    chip1_invert=True)


