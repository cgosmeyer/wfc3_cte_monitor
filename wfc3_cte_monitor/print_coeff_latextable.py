#! /usr/bin/env python

"""
One-off script to help me print coefficient tables for 2016 ISR.
Could be later modified to be more automated.  Retain at least 
as an example.

Authors:

    C.M. Gosmeyer, Nov. 2016

"""

import os
from astropy.io import ascii
from astropy.table import Table
from config import path_to_outputs
from collections import OrderedDict

#-------------------------------------------------------------------------------# 

def print_coeff_latextable(coeff_file_dict):
    """ Prints to terminal the coeffients from text files into a simple 
    latex table format. 
    Prints values to two decimal places.

    Parameters:
        coeff_file_dict : dictionary
            Keys of observation type, values of text filename. 

    """
    coeff_table = Table()
    coeff_table['Mode'] = ['C00', 'C01', 'C02', 'C10', 'C11', 'C12', 'C20', 'C21', 'C22']
    formats = {}
    for colname, coeff_file in zip(coeff_file_dict.keys(), coeff_file_dict.values()):
        data = ascii.read(coeff_file)
        coeff_table[colname] = data['col3']
        formats[colname] = '%0.2e'
        
    ascii.write(coeff_table, format='latex', formats=formats)  


#-------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------# 

if __name__ == '__main__':
    path_to_logflux0 = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf0', 'most_recent', 'cte_vs_logflux')
    path_to_logflux6 = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf6', 'most_recent', 'cte_vs_logflux')
    path_to_logflux12 = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf12', 'most_recent', 'cte_vs_logflux')
    #path_to_logflux18 = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf18', 'most_recent', 'cte_vs_logflux')
    path_to_logflux0_ctecorr = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf0_ctecorr', 'most_recent', 'cte_vs_logflux')
    path_to_logflux6_ctecorr = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf6_ctecorr', 'most_recent', 'cte_vs_logflux')
    path_to_logflux12_ctecorr = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf12_ctecorr', 'most_recent', 'cte_vs_logflux')
    #path_to_logflux18_ctecorr = os.path.join(path_to_outputs, 'finalresults_iraf', 'pf18_ctecorr', 'most_recent', 'cte_vs_logflux')

    # NGC 104
    coeff_files = OrderedDict()
    coeff_files['pf0, short'] = os.path.join(path_to_logflux0, 'cteVSlogflux_ngc104_F502N_s_pf0_ctecorrFalse_r3_mjd57601_coeffs.txt')
    coeff_files['pf0, long'] = os.path.join(path_to_logflux0, 'cteVSlogflux_ngc104_F502N_l_pf0_ctecorrFalse_r3_mjd57601_coeffs.txt')
    coeff_files['pf6, long'] = os.path.join(path_to_logflux6, 'cteVSlogflux_ngc104_F502N_l_pf6_ctecorrFalse_r3_mjd57601_coeffs.txt')
    coeff_files['pf12, short'] = os.path.join(path_to_logflux12, 'cteVSlogflux_ngc104_F502N_s_pf12_ctecorrFalse_r3_mjd57601_coeffs.txt')
    coeff_files['pf12, long'] = os.path.join(path_to_logflux12, 'cteVSlogflux_ngc104_F502N_l_pf12_ctecorrFalse_r3_mjd57601_coeffs.txt')
    #coeff_files['pf18, long'] = os.path.join(path_to_logflux18, 'cteVSlogflux_ngc104_F502N_l_pf18_ctecorrFalse_r3_mjd57601_coeffs.txt')
    
    coeff_files['pf0_ctecorr, short'] = os.path.join(path_to_logflux0_ctecorr, 'cteVSlogflux_ngc104_F502N_s_pf0_ctecorrTrue_r3_mjd57601_coeffs.txt')
    coeff_files['pf0_ctecorr, long'] = os.path.join(path_to_logflux0_ctecorr, 'cteVSlogflux_ngc104_F502N_l_pf0_ctecorrTrue_r3_mjd57601_coeffs.txt')
    coeff_files['pf6_ctecorr, long'] = os.path.join(path_to_logflux6_ctecorr, 'cteVSlogflux_ngc104_F502N_l_pf6_ctecorrTrue_r3_mjd57601_coeffs.txt')
    coeff_files['pf12_ctecorr, short'] = os.path.join(path_to_logflux12_ctecorr, 'cteVSlogflux_ngc104_F502N_s_pf12_ctecorrTrue_r3_mjd57601_coeffs.txt')
    coeff_files['pf12_ctecorr, long'] = os.path.join(path_to_logflux12_ctecorr, 'cteVSlogflux_ngc104_F502N_l_pf12_ctecorrTrue_r3_mjd57601_coeffs.txt')
    #coeff_files['pf18_ctecorr, long'] = os.path.join(path_to_logflux18_ctecorr, 'cteVSlogflux_ngc104_F502N_l_pf18_ctecorrTrue_r3_mjd57601_coeffs.txt')

    print_coeff_latextable(coeff_files)

    # NGC 6791
    coeff_files = OrderedDict()
    coeff_files['pf0, short'] = os.path.join(path_to_logflux0, 'cteVSlogflux_ngc6791_F502N_s_pf0_ctecorrFalse_r3_mjd57598_coeffs.txt')
    coeff_files['pf0, long'] = os.path.join(path_to_logflux0, 'cteVSlogflux_ngc6791_F502N_l_pf0_ctecorrFalse_r3_mjd57598_coeffs.txt')
    coeff_files['pf12, long'] = os.path.join(path_to_logflux12, 'cteVSlogflux_ngc6791_F502N_l_pf12_ctecorrFalse_r3_mjd57598_coeffs.txt')

    coeff_files['pf0_ctecorr, short'] = os.path.join(path_to_logflux0_ctecorr, 'cteVSlogflux_ngc6791_F502N_s_pf0_ctecorrTrue_r3_mjd57598_coeffs.txt')
    coeff_files['pf0_ctecorr, long'] = os.path.join(path_to_logflux0_ctecorr, 'cteVSlogflux_ngc6791_F502N_l_pf0_ctecorrTrue_r3_mjd57598_coeffs.txt')
    coeff_files['pf12_ctecorr, long'] = os.path.join(path_to_logflux12_ctecorr, 'cteVSlogflux_ngc6791_F502N_l_pf12_ctecorrTrue_r3_mjd57598_coeffs.txt')

    print_coeff_latextable(coeff_files)

