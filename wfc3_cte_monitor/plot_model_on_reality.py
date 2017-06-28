"""
Plots the model slopes on the observational slopes vs time.

Authors:

    C.M. Gosmeyer, Aug. 2016

Use:

    >>> python plot_model_on_reality.py

Outputs:

    PNG plots. For each permutation,
    "cteVStime_model_<targname>_<filter>_<exp_length>_pf<flashlvl>_ctecorr<True/False>_r<aperture>.png"
    In "outputs/finalresults/model_on_reality/"

"""

import collections
import glob
import itertools
import numpy as np
import os
import pylab

import config

from astropy.io import ascii

from collections import OrderedDict
from scipy.stats import linregress
from scipy.stats import sigmaclip

from config import program_list
from config import filter_list
from config import target_list
from config import flashlvl_list
from config import fluxbins_lo
from config import fluxbins_hi
from config import path_to_outputs

from set_paths_to_outputs import set_paths_to_outputs


#-------------------------------------------------------------------------------# 

def solve_for_slope(coeffs, mjd, flux, idl):
    """ Calculates slope for set of coefficients, mjd, and flux.
    So it seems that eqs 2-3 in Noeske et al 2012 may have had the time and
    flux parameters swapped: d's should be f's and vice versa.

    Parameters:
        coeffs : list of floats
            List of 9 coefficients.
        mjd : float or int
            Epoch in MJD time.
        flux : float
            Average flux of bin.
        idl : {True, False}
            Set to True if want to plot IDL version's coefficients.

    Returns:
        slope : float
            Empirically derived slope for set of coefficients, mjd, and flux.

    Outputs:
        nothing
    """
    f = np.log10(flux)

    d = float(mjd) - 55400.
 
    """
    # The code from the plotting routine to plot the fits:
    # (This is exactly the same between Python and IDL versions)
    slope = 0
    coeff_count = 0
    for row in range(3):
        for col in range(3):
            print('coeff_count: {}'.format(coeff_count))
            slope += coeffs[coeff_count]*(f**col)*(d**row)
            coeff_count += 1
    """
    if idl:
        # Mixing coefficient numbers because the input file has them ordered 
        # funny and this is the laziest solution.
        slope = coeffs[0] + \
                coeffs[3] * f + \
                coeffs[6] * f**2 + \
                coeffs[1] * d + \
                coeffs[4] * d * f + \
                coeffs[7] * (f**2) * d + \
                coeffs[2] * (d**2) + \
                coeffs[5] * f * (d**2) + \
                coeffs[8] * (d**2) * (f**2)        

    else:
        # Python-derived coefficients.
        # Same code as in the double-for loop in comments above, but 
        # expanded: 
        slope = coeffs[0] + \
                coeffs[1] * f + \
                coeffs[2] * f**2 + \
                coeffs[3] * d + \
                coeffs[4] * d * f + \
                coeffs[5] * (f**2) * d + \
                coeffs[6] * (d**2) + \
                coeffs[7] * f * (d**2) + \
                coeffs[8] * (d**2) * (f**2)


    return slope


#-------------------------------------------------------------------------------# 

def read_coeff_files(targname, filt, exp_length, flashlvl, ctecorr, 
    aperture, basedir='finalresults', idl=False):
    """ Reads in coefficient file of model for the given MJD epoch.
    For Python version's coefficients. 

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
        basedir : string
            Base directory of 'finalresults'. 
        idl : {True, False}
            Set to True if want to plot IDL version's coefficients.
            Much will be hardcoded, because this is a one-off test to 
            show how 'good' the old fits were.
            Also note, there only exists coefficients for the most recent
            epoch (in this case, I think from September 2015)

    Returns:
        coeff_dict : dict of numpy arrays
            Keys of epoch, values of arrays of 9 coefficients.

    Outputs: 
        nothing

    """
    if idl:
        if ctecorr:
            ctedir = '_ctecorr'
        else:
            ctedir = ''

        # change path to .../most_recent_2011 or 2015 to plot the old IDL coefficients from 2012 and 2015 ISRs
        path_to_mostrecent = '/grp/hst/wfc3b/cgosmeyer/Projects/UVIS_CTE_monitor/outputs_idl/finalresults/pf{}{}/most_recent/cte_vs_flux/'.format(flashlvl,ctedir)

        filename = '{}_{}_r{}_coeffs.dat'.format(filt.upper(), exp_length, int(aperture))
        print("coeffs.dat filename: {}".format(filename))

        file_name_search = os.path.join(path_to_mostrecent, filename)

        coeff_dict = collections.OrderedDict()
        if os.path.isfile(file_name_search):
            # Open the file. Read the coefficients.
            coeffs = []
            with open(file_name_search) as f:
                line = f.readlines()
                print(line)
            for s, n in zip(line, range(len(line))):
                if n >= 4:
                    coeffs.append(float((s.split())[-1].split('\n')[0]))

            print("IDL coeffs: ", coeffs)    
            
            for epoch in ['55466', '55629', '55854', '56008', '56131', '56236', '56360',
             '56508', '56673', '56841', '57055', '57221']:
                coeff_dict[epoch] = coeffs

 
    else:
        # Get path.
        path_to_mostrecent = set_paths_to_outputs(
            path_to_outputs, 
            basedir=basedir, 
            flashlvl=flashlvl, 
            ctecorr=ctecorr, 
            mostrecent=True)

        path_to_mostrecent = os.path.join(path_to_mostrecent, 'cte_vs_logflux')

        print('path_to_mostrecent: {}'.format(path_to_mostrecent))

        # Create name of file to search for all available epochs.
        file_name_search = 'cteVSlogflux_{}_{}_{}_pf{}_ctecorr{}_r{}_mjd*_coeffs.txt'\
            .format(targname, filt, exp_length, flashlvl, ctecorr, aperture)   
        print("file_name_search: {}".format(file_name_search))

        files_available = sorted(glob.glob(os.path.join(path_to_mostrecent, file_name_search)))
        print("files available: {}".format(files_available))
     
        # Initialize dictionary where epochs are keys and coeff lists are values.
        coeff_dict = collections.OrderedDict()

        # Loop over all epochs
        for file_name in files_available:

            print('file_name: {}'.format(file_name))

            # Open file if it exists and read coeffs
            if os.path.isfile(file_name):
                data = ascii.read(file_name)
                print(data)
                coeffs = np.array(data['col3']) 
            else:
                print("File does not exist: {}".format(path_to_file)) 

            print("coeffs: {}".format(coeffs))

            # Extract the epoch name.
            epoch = (file_name.split('_mjd')[1]).split('_coeffs.txt')[0]
            print("epoch: {}".format(epoch))

            # Place in dictionary.
            coeff_dict[epoch] = coeffs

    return coeff_dict


#-------------------------------------------------------------------------------# 

def calculate_slopes(coeff_dict, flux, use_latest_coeffs=False, idl=False):
    """ Calculates the slopes using the coeffients in the empirical 
    equation.

    Parameters:
        coeff_dict : dict of numpy arrays
            Keys of epoch, values of arrays of 9 coefficients.
        flux : string or int
            The flux to calculate slope for. Likely the average of the 
            flux bin.    
        use_latest_coeffs : {True, False}
            Plot model using only the latest epoch's coefficients.
            These should be the most precise, since they include the 
            most data.    
        idl : {True, False}
            Set to True if want to plot IDL version's coefficients.

    Returns:
        slopes_model : numpy array
            The slopes calculated with coefficients for each epoch.
        epochs : numpy array
            The epochs corresponding to each slope.

    Outputs:
        nothing
    """

    slopes_model = []
    epochs = []

    coeff_dict_values = list(coeff_dict.values())
    coeff_dict_keys = list(coeff_dict.keys())

    if use_latest_coeffs and coeff_dict_values != []:
        print("----")
        print(coeff_dict_values[-1])
        print(coeff_dict_keys[-1])
        # Use the most recent coefficients for every epoch.
        # Should be the more precise models.
        for coeffs, mjd in zip([coeff_dict_values[-1]]*(len(coeff_dict_keys)), coeff_dict_keys):
            coeffs = np.array(coeffs)
            slope = solve_for_slope(coeffs, mjd, flux, idl=idl)
            print("slope: {}".format(slope))

            slopes_model.append(slope)
            epochs.append(mjd)        

    else:
        # Use the coefficients calculated at each epoch for that epoch.
        # Less precise (less data in the earlier!)
        for coeffs, mjd in zip(coeff_dict_values, coeff_dict_keys):
            slope = solve_for_slope(coeffs, mjd, flux)
            print("slope: {}".format(slope))

            slopes_model.append(slope)
            epochs.append(mjd)

    slopes_model = np.array(slopes_model)
    epochs = np.array(epochs)

    print("slopes_model: {}".format(slopes_model))
    print("epochs: {}".format(epochs))

    return slopes_model, epochs


#-------------------------------------------------------------------------------# 

def read_slope_file(targname, filt, exp_length, flashlvl, ctecorr, aperture, 
    fluxbin, basedir='finalresults'):
    """ Reads in slope file of real slopes over all epochs.

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
        fluxbin : string
            Bin of flux, to be taken average of. 
        basedir : string
            Base directory of 'finalresults'. 

    Returns:
        slopes_cut2 : numpy array
            Slopes of given parameters, including only those with given 
            targname and fluxbin.
        stderrs_cut2 : numpy array
            Slopes' standard errors, including only those with given 
            targname and fluxbin.
        epochs_cut2 : numpy array
            Epochs of the slopes, including only those with given 
            targname and fluxbin.

    Outputs: 
        nothing

    """
    # Get path.
    path_to_mostrecent = set_paths_to_outputs(
        path_to_outputs, 
        basedir=basedir, 
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
        slopes_cut1 = np.array(slopes[np.where(targnames == targname)])
        stderrs_cut1 = np.array(stderrs[np.where(targnames == targname)])
        mjds_cut1 = np.array(mjds[np.where(targnames == targname)])   
        fluxbins_cut1 = np.array(fluxbins[np.where(targnames == targname)])   

        slopes_cut2 = np.array(slopes_cut1[np.where(fluxbins_cut1 == fluxbin)])
        stderrs_cut2 = np.array(stderrs_cut1[np.where(fluxbins_cut1 == fluxbin)])
        mjds_cut2 = np.array(mjds_cut1[np.where(fluxbins_cut1 == fluxbin)])   

        print("slopes: {}".format(slopes_cut2))
        print("stderrs: {}".format(stderrs_cut2))
        print("mjds: {}".format(mjds_cut2))

        return slopes_cut2, stderrs_cut2, mjds_cut2

    else:
        print("File does not exist: {}".format(path_to_file)) 
        return [], [], []


#-------------------------------------------------------------------------------# 

def plot_modelonreality_ctevstime(targname, filt, exp_length, flashlvl, 
    ctecorr, aperture, xlims=[55000,57800], ylims=[-0.1, 1.0], 
    plot_standerrs=True, use_latest_coeffs=False, outloc='', 
    basedir='finalresults', idl=False):
    """ Plots CTE vs time with slopes from the model and from reality.

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
        xlims : list of floats
            The lower and upper bounds of the x-axis (MJD time). 
        ylims : list of floats
            The lower and upper bounds of the y-axis (CTE slope). 
        plot_standerrs : {True, False}
            Set to True if want to plot the standard error as error bars.
            True by default.
        use_latest_coeffs : {True, False}
            Plot model using only the latest epoch's coefficients.
            These should be the most precise, since they include the 
            most data.    
        outloc : string
            Path to the output location. Working directory by default.
        basedir : string
            Base directory of 'finalresults'. 
        idl : {True, False}
            Set to True if want to plot IDL version's coefficients.

    Returns:

    Outputs:

    """

    outname = 'cteVStime_model_{}_{}_{}_pf{}_ctecorr{}_r{}'.format(
        targname, filt, exp_length, flashlvl, ctecorr, aperture)

    fluxbins = ['500-2000', '2000-8000','8000-32000']

    # Can't seem to make dictionaries ordered. So to match the markers
    # and colors with appropriate targnames and filter bins, need make
    # this dictionaries.
    markers = {'observations' : '+', 'model' : 'd'}
    colors = {'500-2000' : 'k', '2000-8000' : 'red', '8000-32000' : 'blue'}

    # Set the figure.
    pylab.figure(figsize=(12.5,8.5))


    for fluxbin in fluxbins:
        print("fluxbin: ", fluxbin)

        # Read in CTE-corrected slopes from observations, to be plotted as-is:
        slopes_observed_ctecorr, stderrs_observed_ctecorr, epochs_observed_ctecorr = \
            read_slope_file(targname, filt, exp_length, flashlvl, 
                ctecorr, aperture, fluxbin, basedir=basedir)

        # Read in non-CTE-corrected slopes from observations, to be plotted as lines:
        slopes_observed, stderrs_observed, epochs_observed = \
            read_slope_file(targname, filt, exp_length, flashlvl, 
                ctecorr=False, aperture=aperture, fluxbin=fluxbin, basedir=basedir)

        # Read in slopes calculated from coefficients.
        # CTEcorr should always be False.
        coeff_dict = read_coeff_files(targname, filt, exp_length, 
            flashlvl, ctecorr=False, aperture=aperture, basedir=basedir, idl=idl)
        #flux = np.mean([float(fluxbin.split('-')[0]), float(fluxbin.split('-')[1])])
        flux = float(fluxbin.split('-')[1])
        print("flux: ", flux)
        slopes_model, epochs_model = calculate_slopes(coeff_dict, flux, 
            use_latest_coeffs=use_latest_coeffs, idl=idl)


        # Plot the observation slopes.
        pylab.plot(sorted(epochs_observed), sorted(slopes_observed), 
            linewidth=3, linestyle='-.', color=colors[fluxbin], alpha=0.6)
        pylab.scatter(epochs_observed_ctecorr, slopes_observed_ctecorr, 
            s=90, marker=markers['observations'], color=colors[fluxbin], 
            label='software {}e-'.format(fluxbin), alpha=0.6)
        if plot_standerrs and stderrs_observed_ctecorr != []:
            pylab.errorbar(epochs_observed_ctecorr, slopes_observed_ctecorr,  
                yerr=stderrs_observed_ctecorr, 
                color='gray', alpha=0.75, marker=None, ls='None')

        # Plot the model slopes.
        pylab.scatter(epochs_model, slopes_model, 
            s=90, marker=markers['model'], color=colors[fluxbin], 
            label='formula {}e-'.format(fluxbin), alpha=0.6)

    print("epochs_model: {}".format(epochs_model))
    print("epochs_observed: {}".format(epochs_observed_ctecorr))
    print("slopes_model: {}".format(slopes_model))
    print("slopes_observed: {}".format(slopes_observed_ctecorr))

    # Complete the plot and save. 
    pylab.xlabel('MJD Date', fontsize=22, weight='bold') 
    pylab.ylabel('CTE loss [flux / 2048 pxl]', fontsize=22, weight='bold')
    pylab.tick_params(axis='both', which='major', labelsize=20)
    # Add lines at major y tickmarks. Must define by hand...
    for ymaj in np.arange(11)*0.1:
        pylab.axhline(y=ymaj,linestyle='--',color='grey')
    title = "{} {} explen'{}' pf{} ctecorr{} ap{}".format(targname, filt, 
        exp_length, flashlvl, ctecorr, aperture)
    pylab.title(title, fontsize=16)
    pylab.xlim(xlims)
    pylab.ylim(ylims)
    pylab.legend(scatterpoints=1, loc='upper left')
    if plot_standerrs and stderrs_observed_ctecorr != []:
        pylab.text(56500, 0.55, 'Standard Error: stddev / sqrt(# sources)')

    if idl:
        pylab.savefig(os.path.join(outloc, '{}_idl.png'.format(outname)), 
            bbox_inches='tight')
    else:
        pylab.savefig(os.path.join(outloc, '{}.png'.format(outname)), 
            bbox_inches='tight')


#-------------------------------------------------------------------------------# 
# The Main.
#-------------------------------------------------------------------------------# 

def plot_model_on_reality_main():
    """ Loops over combinations of targname, filt, exp_length, flashlvl, 
    and ctecorr.

    """
    basedir = 'finalresults_iraf'
    idl = True
    if idl:
        outloc = os.path.join(config.path_to_outputs, basedir, 'model_on_reality_idl')
        flashlvls = [0]
    else:
        outloc = os.path.join(config.path_to_outputs, basedir, 'model_on_reality_python')
        flashlvls = [0, 12] #, 6, 18, 24, 33, 55, 91, 116]

    use_latest_coeffs=True
    aperture = 3
    targnames = config.target_list[:2] 
    filts = ['F502N', 'F606W'] #config.filter_list
    exp_lengths = ['l', 's'] 
    ctecorrs = [False]#[True, False]
    
    for targname in targnames:
        for filt in filts:
            for exp_length in exp_lengths:
                for flashlvl in flashlvls:
                    for ctecorr in ctecorrs:
                        plot_modelonreality_ctevstime(targname=targname, 
                                                      filt=filt, 
                                                      exp_length=exp_length, 
                                                      flashlvl=flashlvl, 
                                                      ctecorr=ctecorr, 
                                                      aperture=aperture,
                                                      use_latest_coeffs=use_latest_coeffs,
                                                      outloc=outloc,
                                                      basedir=basedir,
                                                      ylims=[-0.1, 0.6],
                                                      idl=idl)


# -----------------------------------------------------------------------------


if __name__=='__main__':

    plot_model_on_reality_main()


