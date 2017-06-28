#! /usr/bin/env python

""" 
Initial setup of the uvis_external_cte pipeline.

Sets paths to scripts, outputs, and data. (edits 'config.py')
Need be in the script directory 'detectors/scripts/uvis_external_cte'.
Run when setting up environment for scripts for the first time.
Two steps precede this:
    1. Clone the repo.
    2. (optional) Copy the datasets. By default links to my copies. 
        Here the script can advise which programs are needed. 
See the instructions in the README.

Authors:

    C.M. Gosmeyer, Feb. 2016

Use:
    
    Be in the directory containing the scripts, 'detectors/scripts/uvis_external_cte',
    and do

    >>> python initial_setup.py

Notes:

    Translated from the IDL procedure `initial_setup.pro`.
    Python 2's 'raw_input' is the same as Python 3's 'input' function. 

To Do:
    1. Add creation of database and filling in of Master tables? 
    2. A check that my detectors_tools.photometry package is installed?
    3. Download PA_MAPS from website?
    4. Copy Master catalogs?
"""

from __future__ import print_function

import os


#-------------------------------------------------------------------------------# 

def initial_setup():
    """ Contains initial setup of the uvis_external_cte pipeline.
    """

    print(" ")
    print("Welcome to initial_setup for the WFC3/UVIS External CTE Monitor.")
    print("You will be prompted to enter 4 paths. Afterward, this script will automatically")
    print("set up needed directories and paths. If you need to later edit the paths")
    print("you enter here, edit the program 'set_base_paths.pro' and move directories")
    print("where necessary.")

    # Set path to scripts.
    cwdir = os.getcwd()

    print(" ")
    print("The path_to_scripts is location of uvis_external_cte repo.")
    print("Default path_to_scripts is current working directory:")
    print(cwdir)
    print(" ")
    path_to_scripts = ''
    prompt = '> 1. Type/paste without quotes "." to set the path_to_scripts to CWD. Otherwise type/paste the absolute path: '
    path_to_scripts = input(prompt)

    if path_to_scripts == '.':
        path_to_scripts = cwdir
    else:
        path_to_scripts = path_to_scripts.strip()

    # Set path to data.
    print(" ")
    print("The path_to_data stores all the FLTs and FLCs in proposal ID directories.")
    print("There is no default path.")
    print(" ")
    path_to_data = ''
    prompt = '> 2. Type/paste without quotes an absolute path for path_to_data: '
    path_to_data = input(prompt) 

    path_to_data = os.path.join(path_to_data.strip(), 'data')

    # Set path to outputs.
    print(" ")
    print("The path_to_outputs holds all photometry, tweakreg files, and final plots and coefficients.")
    print("There is no default path. (Nominally make it same as for path_to_data.)")
    print(" ")
    path_to_outputs = ''
    prompt = '> 3. Type/paste without quotes an absolate path for path_to_outputs: '
    path_to_outputs = input(prompt)

    path_to_outputs = os.path.join(path_to_outputs.strip(), 'outputs')

    # Set path to database.
    print(" ")
    print("The path_to_db holds the sqlite database containing master, file info, phot, and results tables for each target.")
    print("There is no default path. (Nominally make it same as for path_to_data.)")
    print(" ")
    path_to_db = ''
    prompt = '> 4. Type/paste without an absolate path for quotes path_to_db: '
    path_to_db = input(prompt)

    path_to_db = os.path.join(path_to_db.strip(), 'database')


    # Edit 'config.py' for these new paths.
    os.rename('config.py', 'config_temp.py') 
    infile = open('config_temp.py', 'r')
    outfile = open('config.py', 'w')

    for line in infile:
        if "path_to_outputs =" in line:
            outfile.write("path_to_outputs = '"+ path_to_outputs + "'\n")
        elif "path_to_data =" in line:
            outfile.write("path_to_data = '" + path_to_data + "'\n")
        elif "path_to_scripts =" in line:
            if path_to_scripts != ".":
                outfile.write("path_to_scripts = '" + path_to_scripts + "'\n")
            else:
                outfile.write(line)
        elif "path_to_db =" in line:
            outfile.write("path_to_db = '"+ path_to_db + "'\n")
        else:
            outfile.write(line)

    infile.close()
    outfile.close()
    os.remove('config_temp.py')

    def create_dir(basepath, dirname):
        fullpath = os.path.join(basepath, dirname)

        if not os.path.isdir(fullpath):
            os.mkdir(fullpath)
            print('Creating {}'.format(fullpath))

        return fullpath

    # Create the needed '/data' base directories. 
    create_dir(path_to_data, '')
    create_dir(path_to_data, 'newdata')
    create_dir(path_to_data, 'PA_MAPS')
    create_dir(path_to_data, 'master_cat')

    # Create the 'database' directory.
    create_dir(path_to_db, '')

    # Create the needed '/outputs' base directories. 
    create_dir(path_to_outputs, '')
    path_to_crmasks = create_dir(path_to_outputs, 'crmasks')
    path_to_finalresults = create_dir(path_to_outputs, 'finalresults')
    path_to_photom = create_dir(path_to_outputs, 'photom')
    path_to_tempfiles = create_dir(path_to_outputs, 'tempfiles')
    path_to_tweakreg = create_dir(path_to_outputs, 'tweakreg')


    # Create a README.txt for each '/outputs' sub-directory.

    readme1 = open(os.path.join(path_to_crmasks, 'README.txt'), 'w')
    readme1.write("The cosmic-ray masks for each image, for the first and second SCI extensions.\n")
    readme1.write("Generated using AstroDrizzle in 'run_adriz.py.'")
    readme1.close()

    readme2 = open(os.path.join(path_to_finalresults, 'README.txt'), 'w')
    readme2.write("The final plots and coefficient files.")
    readme2.write("Directories are time-stamped. The most recent of all versions is copied to 'most_recent'.")
    readme2.close()

    readme3 = open(os.path.join(path_to_photom, 'README.txt'), 'w')
    readme3.write("The *coo and *mag files, and any files related to the photometry.")
    readme3.close()

    readme4 = open(os.path.join(path_to_tempfiles, 'README.txt'), 'w')
    readme4.write("A space to temporarily generate files before they are moved or deleted.")
    readme4.close()

    readme5 = open(os.path.join(path_to_tweakreg, 'README.txt'), 'w')
    readme5.write("The tweakreg output files to be used to identify sources and validate")
    readme5.write("tweakreg solutions.  The files are:\n")
    readme5.write("* ipppssoot_xy_catalog.match\n")
    readme5.write("* ipppssoot_flc_catalog_fit.match\n")
    readme5.write("* ipppssoot_flc_tweakreg.log\n")
    readme5.write("* ipppssoot_hlet.fits\n")
    readme5.write("* hist2d_ipppssoot_flc.png\n")
    readme5.write("* residuals_ipppssoot_flc.png\n")
    readme5.write("* vector_ipppssoot_flc.png")
    readme5.close()


    print(" ")
    print("Initial setup of UVIS External CTE Monitor complete.")
    print(" ")
    print("Before you can run the pipeline, you still need do the following: ")
    print("1. Download the Pixel Area Maps from the website (see README).")
    print("2. Copy the Master catalogs from cgosmeyer (see README).")
    print("3. Create the database and populate the Master tables (see README).")
    print("4. Install the 'detector_tools' repo, primarily for 'photometry' subpackage.")

#-------------------------------------------------------------------------------# 
# The Main.
#-------------------------------------------------------------------------------# 

if __name__=='__main__':
    initial_setup()