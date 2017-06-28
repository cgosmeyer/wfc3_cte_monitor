WFC3/UVIS External CTE Monitor
==============================

Author: C.M. Gosmeyer

Created: 10 Apr. 2015

Last Update: 12 Feb. 2017

**NOTE:**
**Originally located in private repository `spacetelescope/detectors/scripts` as `uvis_external_cte`.**
**Written for the Hubble WFC3 team at Space Telescope Science Institute.**
**Some modificiations have been made to this version to keep file structures hidden.**
**All Hubble data used in this monitor are public.**


-------------------------------------------------------------------------------


Contents
++++++++

0. Introduction

1. Instrument Science Reports and Relevant Documents

2. Proposals

3. Running the Pipeline

   3.1 Initial Setup
   3.2 Procedure for Ingesting ALL the Data
   3.3 Procedure for General Run of the Pipeline
   3.4 Procedure for Processing the 180-Degree ('Same-Chip') Dataset
   3.5 Procedure for Switching Between Python and IRAF Photometry
   3.6 Example Use Cases of the Pipeline

4. Care and Feeding of the Database

   4.1 Procedure for Adding/Renaming a Column
   4.2 Procedure for Removing a Column
   4.3 Procedure for Adding a Table
   4.4 Procedure for Re-Setting and Re-Populating Database
   4.5 Help! The Database Locked Me Out!

5. Summary of Scripts and Modules

   5.1 Pipeline
   5.2 Database
   5.3 For Tests and ISRs

6. Issues and Future Improvements


-------------------------------------------------------------------------------


0. Introduction
+++++++++++++++

This document is meant as a guide for setting up and running the pipeline 
of the external Charge Transfer Efficiency (CTE) monitor. All scripts are 
descended from Kai Noeske’s “PHOTCODE” suite. The originals can be found 
by going way back to the first commits in Quicklook's
`spacetelescope/detectors/scripts/uvis_external_cte`. 

This program monitors the CTE using observations of star clusters (hence 
the 'external' bit), and via these data, produces a table of coefficients 
that can be used to correct for CTE loss given an observation’s epoch, 
distance from amplifiers, source flux, and exposure length. 

This program has many, many parameters to keep in mind.
* Two star clusters of different densities are observed: sparse NGC 6791 
and dense NGC 104 (47 Tuc).
* Two filters were used until Cycle 20: F502N and F606W. 
* Following Cycle 20 (inclusive) all observations are taken only in F502N 
so there will be time to take post flash images of each field. In 
addition, it was deemed unnecessary to continue F606W because post flash 
could provide, as the two filters did, a variety of backgrounds. Non-post 
flash images are still taken for continuity. 
* Two flash levels are sampled: 6 and 12 e-/pxl, of both clusters.
* For dense NGC104, flash exposures 1.24, 1.59, 13, and 21.5 are sampled 
to get a range of background levels. (Filter 502N has effectively no 
background.) This gives us effectively flash levels of 33, 55, 91, and 
116 e-/pxl.
Starting in 14378, we sample finer gradients of flash levels: 18 and 24 
e-/pxl.
* Flash 12 e-/pxl is most used in the analysis, as it is the one recommended 
to observers for general use.
* All images have the pixel-based CTE correction run over them. Every FLT, 
then, should have a corresponding FLC. 
* Photometry is run on apertures 3, 5, and 10 pixels. Aperture 3 is the 
one used in most of the analysis use because it is most reliable - any 
larger aperture allows in too many contaminating sources.

See ISRs and Proposals below for further background information.  Again,
this document's focus is on how to run the pipeline, not on the monitor's 
observation setup or on analysis of results. 


-------------------------------------------------------------------------------


1. Instrument Science Reports and Relevant Documents
++++++++++++++++++++++++++++++++++++++++++++++++++++


"The Efficacy of Post- Flashing for Mitigating CTE-Losses in WFC3/UVIS Images."
J. Anderson, J. MacKenty, S. Baggett, and K. Noeske, K., 2012
http://www.stsci.edu/hst/wfc3/ins_performance/CTE/ANDERSON_UVIS_POSTFLASH_EFFICACY.pdf

WFC3 ISR 2012-09: "WFC3 UVIS Charge Transfer Efficiency October 2009 to October 2011"
Noeske et al., 08 Jun. 2012
http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2012-09.pdf

ACS ISR 2012-05: "A new accurate CTE photometric correction formula for ACS/WFC"
Chiaberge, M., Oct. 2012
http://www.stsci.edu/hst/acs/documents/isrs/isr1205.pdf

WFC3 ISR 2015-03: "WFC3/UVIS Charge Transfer Efficiency 2009 - 2015"
S. Baggett, C. Gosmeyer, K. Noeske, 31 March 2015
http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2015-03.pdf

WFC3 ISR 2016-17:  "WFC3/UVIS External CTE Monitor: Single-Chip CTE Measurements"
C. Gosmeyer and S. Baggett, 14 Dec. 2016
http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2016-17.pdf

WFC3 ISR 2017-09: "WFC3/UVIS External CTE Monitor: 2016 Updates on Coefficients and Analysis Pipeline"
C. Gosmeyer and S. Baggett, 2017
http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2017-09.pdf

WFC3 Instrument Handbook for Cycle 23, "6.9 Charge Transfer Efficiency"
http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c06_uvis10.html


-------------------------------------------------------------------------------


2. Proposals
++++++++++++

(Based on Table 1 of ISR 2015-03)

Program	Cycle	Obs Dates			Filters		Targets
-------------------------------------------------------------------------------
11924	17	Oct 2009,Mar 2010,Sep 2010	F502N, F606W	NGC6791
http://www.stsci.edu/hst/phase2-public/11924.pdf

12348	18	Sep 2010			F502N		NGC104
http://www.stsci.edu/hst/phase2-public/12348.pdf

12379	18	Nov 2010,Mar 2011,Apr 2011	F502N, F606W	NGC6791, NGC104
http://www.stsci.edu/hst/phase2-public/12379.pdf

12692	19	Oct 2011,Mar 2012,Jun/Jul 2012	F502N, F606W	NGC6791, NGC104
http://www.stsci.edu/hst/phase2-public/12692.pdf

13083	20	Nov 2012,Mar 2013,Jul/Aug 2013	F502N		NGC6791, NGC104
http://www.stsci.edu/hst/phase2-public/13083.pdf

13566	21	Jan 2014,Jul 2014		F502N		NGC6791, NGC104
http://www.stsci.edu/hst/phase2-public/13566.pdf

14012	22	Feb 2015,Jun 2015		F502N		NGC6791, NGC104
http://www.stsci.edu/hst/phase2-public/14012.pdf

14378 23  Jan 2016, Jul-Aug 2016  F502N,  NGC6791, NGC104
http://www.stsci.edu/hst/phase2-public/14378.pdf

14541 24 Jan 2017, ...
http://www.stsci.edu/hst/phase2-public/14541.pdf

Notes:

* 12348 was not a monitor proposal; its purpose was to test the charge injection
  (CI) mode as a means of mitigating CTE. But we use some of the images
  because one of our standard clusters, NGC104, was observed in non-CI mode. 
  (The code is capable of sorting through CI and non-CI modes. For now it just
  tosses the CI mode images out, but this could always be modified in 
  run_results.pro if you wish.)

* 12379, Visit 13 not only was imaged with CI but also had a guide star failure. 
  (See http://www.stsci.edu/cgi-bin/get-visit-status
id=12379&markupFormat=html&observatory=HST)

* 12692, Visits 10 and 11, consists of the 180-degree dataset, observed with
  a non-nominal field in ngc6583.

* 14012, Visit 05 is not part of the normal monitor (the data are binned for a 
  separate experiment) and should be removed from the proposal directory.


-------------------------------------------------------------------------------


3. Running the Pipeline
+++++++++++++++++++++++

The scripts are written in Python.  They were originally in IDL, Fortran, 
and IRAF, written by V. Platais (geometric distortions) and K. Noeske 
(everything else).  They were refurbished by C.M. Gosmeyer in Spring-Summer 
2016.  

You might find the originals by going back in the master branch's history.
The new pipeline was merged on 30 Nov. 2016.

The pipeline is capable of performing on many combinations of data modes 
(postflash, pixel-based cte-corrected, aperture radius sizes, etc.). In the
subsections below we describe the most general cases for running the pipeline.
We will not go into the inner workings of each module and function here. 
Please refer to their doc strings and comment lines. 


3.1 Initial Setup
-----------------

A new user of the scripts should first complete all steps in this section. 

I have checked that the pipeline (with Python-based photometry) works in 
Astroconda with Python 3.  You should try setting up in a Python 3 environment
initially.
To run the pipeline with IRAF photometry, you'll need a Python 2.7 environment
and you may need to wrestle with it to get IRAF running... 

1. Clone the `wfc3_cte_monitor` repo, to wherever convenient. The scripts
   are independent of the input data and the output file’s locations.
   Run the `setup.py`:

   >>> python setup.py develop

2. In your `wfc3_cte_monitor\wfc3_cte_monitor` directory run the script

   >>> python initial_setup.py

   and follow the prompts. You should now have a “data” and an “outputs” 
   directory where you specified.

   In addition, your `config.py` will be modified to contain the paths you 
   specified to the script. Have a look at `config.py` to check all paths are 
   correct. If you ever change the locations of your repo, "data", or 
   "outputs" directories, you only need modify this script.

   The “data” directory will eventually contain *flt.fits and *flc.fits 
   files sorted into directories by proposal number. Each Cycle a new 
   directory will be added.  These proposal directories will be generated 
   when new data is processed with `run_tweakreg.py`.

   "data" should contain the following directories (initially empty).
       master_cat - contains the Master catalogs
       newdata - the holding directory for new data
       PA_MAPS - contains the pixel area maps
   
   The “outputs” directory consists of the following four sub-directories.
   They initially contain only READMEs describing their contents.
       finalresults - contains final plots and coefficients.
       photom - contains coo, mag, and diagonistic plots from photometry.
       tempfiles - holding folder for temporary files.
           (this may be degenerate, as nothing in the pipeline currently 
           needs to store temp files)
       tweakreg - contains intermediate files and diagonistic plots from
           TWEAKREG. 

3. Clone cgosmeyer's `photutils_plus` repo. Go to its top level and install
   as

   >>> python setup.py develop

   You should then be able to import the 'photutils_plus' package.  Check this 
   in ipython.

   python> from photutils_plus import photutils_plus.

4. Likewise, clone cgosmeyer's `analysis_tools` repo and run its `setup.py`.

5. Obtain the Pixel Area Maps.
    
   Into "data/PA_MAPS" download the pixel area maps for UVIS chips 1 and 2 
   from http://www.stsci.edu/hst/wfc3/pam/pixel_area_maps

   You should now have 

       UVIS1wfc3_map.fits
       UVIS2wfc3_map.fits     

6. Copy the Master catalogs for the three clusters.  Currently the *cat files
   and the drizzled images used to create them are located on C.M. Gosmeyer's
   central store.  

   In "data/master_cat" copy directory from [C.M. Gosmeyer will provide you with path]

   The pipeline will expect to find the *cat files in "data/master_cat/final".

7. Set up the SQLITE database.  

   If you are in a hurry and want an already filled-in sqlite database file
   "uvis_external_cte.db", copy it from C.M. Gosmeyer's central store to 
   your "database" directory. 
   Note that if you do this, you may not be able to repopulate the FileInfo
   and Phot tables with the `database_reset.py` functions because they rely
   on the paths to the photometry files stored in the database. You could 
   probably get around this with some hacks.  In any case, if you are 
   taking over this monitor, you probably should populate the database from
   scratch in your own environment at some point.

   To create the "uvis_external_cte.db" file and initialize the tables in
   the location you specified in `config.py` run 

   >>> python database_interface.py   

   To fill in the Master tables run

   >>> populate_master_tables.py

   Check that the database is correctly filled in. Go to the location of
   the SQLITE database and enter it with the following command.

   >>> sqlite3 uvis_external_cte.db

   Then do

   sqlite> .tables

   You should see all the following tables.

    ngc104_fileinfo   ngc104_results    ngc6583_phot      ngc6791_master  
    ngc104_master     ngc6583_fileinfo  ngc6583_results   ngc6791_phot    
    ngc104_phot       ngc6583_master    ngc6791_fileinfo  ngc6791_results 

   Then to check that each table has all the correct columns, do, for 
   example,
    
   sqlite> pragma table_info(ngc104_master);

   You should see all the columns that were defined in the files located 
   in "detectors/scripts/uvis_external_cte/table_definitions"

   Finally check that the Master tables are filled in. For example,

   sqlite> select id,xpix,ypix,ra,dec from ngc6583_master;

   You should see columns printed containing the catalog information
   from the .cat files located in "master_cat/final/".

   To quit SQLITE, do

   sqlite> .exit

   Note that there will be *no* entries in the columns until you ingest
   the data, as outlined in Section 3.2. 

8. You should now to ready to run the scripts! Follow the procedure below
   for ingeseting ALL the data if you haven't yet copied the data into
   "outputs".

* A note about the `config.py`*  
You should *not* commit this once you add your personal paths.  
In this file you will also find static lists of flux bin ranges, proposal
IDs, target names, filters, and flash levels.  If you ever need these lists
just import them from here so you know they will always be consistently the 
same across all functions (and hopefully all correct!).
An example import of the targets would be 

    from config import target_list

Note that all the looping lists and switches in the pipeline are either set 
in the config.py or are arguments (some of them with default values) in the
primary wrapper, `run_uvis_external_cte.py`. 
The idea is that, to change what subsets of data it is legal for the 
pipeline to chew on, modify the config.py. To constrain the pipeline from
looping over and doing EVERYTHING, use the command-line switches parsed by
`run_uvis_external_cte.py`. 
This idea has been imperfectly implemented as yet; see Section 6. 


3.2 Procedure for Ingesting ALL the Data
----------------------------------------

1. If you do not already have them in your "data" directory, retrieve 
   all the FLCs and FLTs from MAST for all proposal numbers listed above
   in Section 2 and store them in "data/newdata".
   Alternatively, copy all the RAWs from Quicklook to "data/newdata" 
   and run calwf3 over them. You want all the images to be calibrated 
   with the same version of the calibration pipeline.

2. Procede as in Step 3 of "3.3 Procedure for General Run of the Pipeline".
   All images should be copied into subdirectories named for the proposal
   numbers in the "data" directory.

3. If you are running all data for the first time (i.e., you have not yet 
   filled in the FileInfo or Phot tables), do

   >>> python run_uvis_external_cte.py --pr 'all' 

   This step also works if you want to re-fill in the database, i.e.,
   if you found an error in photometry. All entries will get overwritten.

   The script `run_image_extraction.py` is smart enough to skip over bad
   or irrelevant datasets within the proposals, such as those described
   above in Section 2. See the script's functions' docstrings for more 
   details.

4. Wait. Filling in the database could take dayyys unless you can figure 
   out how to multi-process things... The slow-down seems to be writing
   the photometry results to the database.
   See Section 6 for more details.

5. Make a copy of your filled database in a seperate, safe directory (even 
   better, on a different partition or machine) so you can fall back
   to it should anything go horribly wrong with the active database.  
   Back it up as needed.

6. Read Section 3.3 for more information on individual steps.


3.3 Procedure for General Run of the Pipeline
---------------------------------------------

Follow this procedure if you have a new epoch to ingest.
If starting a new Cycle, be sure that the proposal list in `config.py` 
is updated. (And, of course, update this file's "Programs" section!)

1. Copy the FLTs and FLCs from MAST or Quicklook directories to 
   "data/newdata".
   If you are copying ALL the data to be run through scripts for the first
   time, you should either go through MAST or run calwf3 on the RAWs yourself 
   so you can ensure all files will be calibrated with the same and most 
   recent version of calwf3.

2. Do a quick quality check of all the data. (Open in DS9 or look at
   the Quicklook JPEGs.) Make sure all expected data for a given visit 
   is there. You can check this by comparing to the proposal pdf.
   You might also be clever and write a script that checks all expected
   data is in fact there.

3. In "data/newdata" run 

   >>> python run_tweakreg.py

   This will correct the WCS of each image so that the photometry scripts 
   can match the stars exactly to the Master catalog using RAs and Decs.
   Don't be discouraged if this is a bit slow!

   Once the script is complete, go to the directory "outputs/tweakreg"
   and look through the diagonistic plots.
   * hist2d_*_flc.png - 2d histogram bins should be centered.
   * residuals_*_flc.png - scatter points should be evenly distributed
     around line without systematics.
   * vector_*_flc.png - vectors should be random and evenly distributed.

   The script will move the FITS files to the appropriate proposal directory
   in "outputs/data". Now you should be ready to run the pipeline.

4. Run the pipeline wrapper `run_uvis_external_cte.py`. It does the 
   following over CTE-corrected on non-CTE-corrected images of NGC104 and 
   NGC6791 for all available flash levels and exposure lengths (long and 
   short).
   * Matches sources to the Master catalog.
   * Performs photometry on the matched sources.
   * Matches the image pairs.
   * Measures CTE decline per image pair by slopes of Chip 1/Chip 2 flux-
     ratios vs Chip 2 y-positions.
   * Creates plots of CTE decline vs log10 flux, fitting polynomials to each 
     epoch, obtaining coefficients for the aperture photometry correction 
     formula. Displays these plots on the Quicklook website.
   * Creates plots of CTE decline vs MJD. Displays these plots on the 
     Quicklook website.
   * Fills in the FileInfo, Phot, and Results tables in the database. 

   The wrapper's parameters are by default set to assume the most recent 
   ('last') proposal is being ingested. So to run the pipeline in this case, 
   you only need specify one parameter,

   >>> python run_uvis_external_cte.py --pr 'last' 

   Read the doc strings of `run_uvis_external_cte.py` for more options.
   You can also see all the parameter options on the commandline by

   >>> python run_uvis_external_cte.py -h

   The pipeline can be run over a specific proposal or over all the 
   proposals, as well as many other variations of plots created and options 
   of whether to insert new data into the database, whether to copy new 
   files to "automated_outputs" for display on Quicklook website, etc.

   Note that the 'last' option will re-run *all* the visits through the
   pipeline for that most recent proposal, which you may not want if you
   are time-constrained and you already ran the first visits through.
   It's not bad if a visit gets re-run.  It will just take extra time. 

   This functionality could be changed. Really is up to you; it only exists 
   as it is so you don't necessarily need to remember the proposal name to 
   run the pipeline on the newest data nor remember whether you already ran 
   it with the previous visit in that proposal.  It was made to be a dummy-
   proof. There may have been instances in the old pipeline days where the
   runner pushed the most recent data through and could not figure out why 
   that proposal's first visit from 6 months ago was missing from the plots 
   (she forgot she hadn't run that one through yet...) or why the pipeline 
   wasn't generating photometry files for the most recent proposal (she had 
   been specifying the name of an older proposal for the pipeline to run...)
   
   Anyway! You can instead specify the exact proposal number and exact 
   visit number(s) like in this example:

   >>> python run_uvis_external_cte.py --pr '14012' --v 3 4


3.4 Procedure for Processing the 180-Degree ('Same-Chip') Dataset
-----------------------------------------------------------------

The 180-degree dataset is from proposal 12692, visits 10 and 11,
imaging the non-nominal star cluster NGC 6583. Its purpose was
to provide a dataset from which the CTE loss of Chip 1 and Chip 2
individually could be measured. See the proposal and ISR 2016-17
for more information.

Because it uses its own star cluster, this dataset has its own 
tables in the database: ngc6583_master, etc. 

Run the FLTs and FLCs through the tweakreg scripts as you normally would 
for the other datasets. Then to run the pipeline on these data, simply do

>>> python run_uvis_external_cte.py --pr '180'

Outputs will be ingested into the database and the files will be
sorted just like the other datasets. For example, outputs/finalresults/pf0/most_recent/ngc6583_F502N_l

Note that there will only be flashlvl=0, filter=F502N (long and short 
exposures), and filter=F606W (long exposure only).

Included within the `uvis_external_cte_plots.py` module is a special
180-degree dataset plot from the function `plot_180_slope_vs_expt_setup`,
which is run automatically from the pipeline in `run_outputs.py` when 
you set from the commandline, --pr '180'. This function creates plots in 
outputs/finalresults/pf0/most_recent/180cte_expt/

In addition to the nominal plots (which by themselves really aren't
too enlightening) there is a special plotting script for this dataset:

* cte180test_plots.py -- Overplots the 180-degree dataset onto the 
    nominal NGC 104 data from the same 12692 proposal onto a CTE loss
    vs log10 flux plot. This was the primary plot of ISR 2016-17.

A few words about how the 180-degree dataset is run through the pipeline.
You can see that throughout there are '180' if-else switches. But
really the photometry, source-to-master-catalog matching, etc are all
the same.  The one spot where there's a substantial difference is in how 
the image pairs are matched.  

Right now the image-pair matching is a bit hokey. Go to `run_image_extraction.py`
in the function `create_param_dict`.  There is the line

  if '6583' in targname:

This procedes to assign images to chip numbers based on hard-coded strings, 
which are unique letters in each of the image's filenames.  As a one-off test, 
it's not the worst thing ever...  But if the monitor starts taking these 
images regularly, you should come up with something more automated!

The matched pairs are later retrieved using the database querying function
`query_for_180pair` from module `database_queries.py`. You can see it 
utilized in the 180-degree plotting functions described above.
It relies on the chip numbers assigned in `run_image_extraction.py`.
Matching pairs uniquely is tricky in the 180-degree test. The nominal
dataset's pairs can be matched based solely on the different chip 
numbers (where all other parameters match). But in the 180-degree
test pairs *have the same chip number*.  Thus the `query_for_180pair`
function relies on the observation setup to provide pairs whose
filenames come in alphabetical order, should all other parameters
be the same - for example, there are two pairs of long exposures 
that I wanted to ensure would be unique. Except for consistency, there 
was no other reason to match one image to the other. (Relying on 
alphabetical order was the lazy way to do this. Could do a more 
sophisticated retrieval of the 'obstime' from the headers or whatever.)



3.5 Procedure for Switching Between Python and IRAF Photometry
--------------------------------------------------------------

I actually have two SQLITE databases which are identical except that
one contains photometry from Python's `photutils.aperture_photometry`
function and the other from IRAF (pyraf-wrapped) `phot` procedure.

To switch between them, you need do two things.  

First, you need change the name of the database in the `config.py` file.

For my Python-based database:
    db_name = 'uvis_external_cte.db'

For my IRAF-based databse:
    db_name = 'uvis_external_cte_iraf.db'

Second, when you run the full pipeline (with the intent to do photometry), 
you need switch IRAF on with the following keyword in order for IRAF-readable 
coordinate files to be generated, photometry be performed with the `iraf.phot` 
function, and for outputs to be placed in the proper *_iraf subdirectories 
of "outputs".

>>> run_uvis_external_cte.py --irafphot 'y'

By default this option is 'n'.

If you are just making plots, you don't need to specify IRAF on the 
commandline because it will extract data from the active database,
which if you modified `config.py` correctly, should be the IRAF one.
(One could make this more dummy-proof, but works for now...)

See Section 6 for an explanation why I stuck to IRAF photometry for the
time being.


3.6 Example Use Cases of the Pipeline
-------------------------------------

Below we give examples for running the scrips in several scenarios.

* If re-running photometry and databse ingestion over just a specific 
  program (ie, 12379) for aperture 3 and only want to recalculate the
  slope plots,

  >>> python run_uvis_external_cte.py --pr '12379' --ap 3 --pl 'ratio_ypos'


* If want to do all photometry and plots for 180-degree dataset,

  >>> python run_uvis_external_cte.py --pr '180' --pl 'ratio_ypos' '180cte_expt' --ao 'n'

And so on.   

Also, this pipeline carries with it an extensive database of stars uniquely
identified over many observing modes through over six years of continuous 
observation. Many side projects could result from this!  



-------------------------------------------------------------------------------

4. Care and Feeding of the Database
+++++++++++++++++++++++++++++++++++

The database is SQLITE, basically a file with the extension DB.  It can be
copied, moved, or whatever. Each database 'file' is made up of tables (e.g., 
master_ngc104, etc.). Each table is made up of columns (e.g., exptime, filter, 
etc.). Entries can be overwritten. 

The columns for each table are defined in the text files contained in the 
subdirectory "table_definitions". 

See Section 5 for a list of all modules and scripts relevent for interacting
with the database. 


4.1 Procedure for Adding/Renaming a Column
------------------------------------------

1. In the "table_definitions" subdirectory of this repo, open the table
   for which you want to add a column. Then type the name of the column
   followed by the datatype (String, Float, or Integer).

2. In `database_update.py`, go to the appropriate `results*dict` function
   and add the lists, etc., which all should be obvious in context.
   The hardest one to modify is `return_phot_dict_list`, since that has 
   so many columsn already for the different photometry apertures and
   because you need enter them differently for IRAF vs Python-based 
   photometry.

3. Copy the database for safety.
   >>> cp uvis_external_cte.db uvis_external_cte_copy.db

4. Enter SQL to drop the table for which you want to add/rename a column.
   This is permanent!  Be sure you saved a copy of the database!!

   >>> sqlite3 uvis_external_cte.db

   sqlite> drop table table_name;

   sqlite> .exit

5. Re-generate the table. It should now contain the new/renamed column.

   >>> python database_interface.py  

   Check that the table and its new/renamed column is there by entering 
   SQL again.

   >>> sqlite3 uvis_external_cte.db

   sqlite> .tables

   sqlite> pragma table_info(table_name);

   sqlite> .exit

6. Re-populate the table.  Feel free to modify this step as needed.  It
   will depend on your column. You might be able to just use the 
   `database_reset.py` script to repopulate columns and write another 
   small script to populate the new column.  Or you may have to re-run
   the pipeline entirely, especially if your column changes have to do
   with photometry.

7. If you are satisfied that nothing went horribly wrong, delete the
   copy of the database (presumably you'll have another backup
   somewhere in a different safe directory anyway!).


4.2 Procedure for Removing a Column
-----------------------------------

1. In the "table_definitions" subdirectory of this repo, open the table
   that contains the column and delete its entry.

2. In `database_update.py`, remove all references to it in the appropriate 
   `results*dict` function.

3. Procede as in Section 4.1 to drop the table and repopulate it without
   the column.


4.3 Procedure for Adding a Table
--------------------------------

Likely any new table you create will need have counterparts for each
target. For examle: master_ngc104, master_ngc6791, master_ngc6583.
Each instance of a table class will have the same column names and use 
the same updating and retrieving functions.

1. Create a new file in the "table_definitions" subdirectory of this repo.
   The name of the file must be the base name of the new table. Look into
   the existing files for how to name the columns within the file.

2. In `database_interface.py`, create an "orm" function for the table class.
   Then intilize it for each target, following the naming convention you'll
   see, under the ``# Initialize classes`` line.
   Finally, add the new table instance names to the ``__all__`` list at the
   top of the module.

3. In `database_update.py` create an "update" and "return" function for
   the table class.
   Add imports from `database_interface` for your new table's instances 
   to the top of the file. For example,

       from database_interface import NGC104_<newtable>
       from database_interface import NGC6791_<newtable>

   Note that in any of the modules that need access to the tables in
   order to query or update them, you'll need add these imports.

4. Run `database_interface.py` to generate the new tables.

   >>> python database_interface.py  

   Go into SQL to check that the new tables are there with the correct 
   columns.

   >>> sqlite3 uvis_external_cte.db

   sqlite> .tables

   sqlite> pragma table_info(table_name);

   sqlite> .exit

5. You may consider adding a "reset" function to `database_reset.py`.  
   May or may not be viable for your new table class.


4.4 Procedure for Re-Setting and Re-Populating Database
-------------------------------------------------------

You *could* drop all the tables as described in step 4 of Section 4.1. 
But it's probably easiest just to delete (in reality, stash away in case
the worst happens) the *db file and re-generate the database as described
in Step 5 of Section 3.1.  Then you can, if you wish, re-populate the
Master, FileInfo, and Phot tables using the data from the old database.
Go into `database_reset.py` and modify the Main in order to reset the
tables you desire.

>>> python database_reset.py

Note that there is not currently a function for reseting the Results tables.
You'll still need to run the pipeline like this:

>>> python run_uvis_external_cte.py --pr 'all' --dofind 'n' --dophot 'n' --pl 'ratio_ypos'

The Master tables always need be generated first.  If starting from scratch
and you don't want to re-populate with old data, just repopulate the Master
tables with:

>>> python populate_master_tables.py


4.5 Help! The Database Locked Me Out!
-------------------------------------

Sometimes the database 'locks' when it the pipeline crashes mid-injestion.
Fortunately unlocking it is a two-command procedure. Go to where the database
is located and do

>>> fuser uvis_external_cte_iraf.db

You will get a message about what process locked the database. For example,

    uvis_external_cte_iraf.db: 1234

To unlock the database, you just need to kill the process (here, '1234'):

>>> kill -9 1234



-------------------------------------------------------------------------------


5. Summary of Scripts and Modules
+++++++++++++++++++++++++++++++++

5.1 Pipeline 
------------

* config.py -- filled by initial_setup.py. Used for setting global lists
   and paths, for import.

* initial_setup.py -- run this to initialize directories for pipeline.

* run_adriz.py -- runs AstroDrizzle to create CR mask and to obtain a 
    background measurement in MDRIZSKY keyword. It will be up to results
    of your testing whether to officially make this the second prep step
    before running actual pipeline.

* run_image_extraction.py -- performs the setup of the images. This includes 
   pulling out all an storing all relevant header information from the images 
   and performing source-finding, catalog-matching, and photometry on them. 
   From here the FileInfo and Phot tables get filled. 

* run_outputs.py -- measures CTE and outputs various plots used for analysis.
   From here the Results table (CTE measurements) gets filled.

* run_tweakreg.py -- runs TweakReg to correct WCS in new image's headers
   to WCS in master image of that field. The first prep step before 
   running pipeline itself.

* run_uvis_external_cte.py -- the pipeline-running script. 

* set_paths_to_outputs.py -- handles all logic to name paths to outputs 
   in which we want flashlvl-ctecorrected-timestamp ordered directories. 

* uvis_external_cte_plots.py -- contains plot-creation functions. 


5.2 Database
------------

* database_interface -- Module to interface and connect to the database.  
  When run as script on command line, creates initial empty database.

* database_queries -- Module containing the database queries.

* database_reset -- Module for re-populating the database's FileInfo, 
  Master, and Phot tables, without having to redo the photometry calculations.
  Can be run as a script by modifying the Mean.

* database_update -- Module for interfacing with database to update the 
  various tables either by inserting new records or updating existing ones.

* populate_master_tables -- Populates all the Master tables.


5.3 For Tests and ISRs
----------------------

* compare_idl_python.py -- overplots IDL pipeline CTE measurements on
   Python pipeline measurements. For ISR 2017-xx.

* cte180test_plots.py -- overplots 180-degree log10 flux vs CTE degredation
   on nominal data from same proposal. For ISR 2016-17.

* photom_tests_modeldata.py -- for comparing iraf.phot to photutils.
    aperture_photometry on a simulated 2-d gaussian star.

* plot_model_on_reality.py -- plots the model CTE-measuring slopes on the 
    observational slopes vs time. For ISR 2017-xx.

* print_coeff_latextable.py -- prints a base LaTeX table from coefficient
    text files, for ease of incorporating into ISR 2017-xx.


-------------------------------------------------------------------------------


6. Issues and Future Improvements
+++++++++++++++++++++++++++++++++

These can be used as exercises for the new monitor lead to get used to the scripts
and database.  They are in no particular order of importance. 

1. Master catalogs could be better. 
   There's an issue with the wide-band F606W data. The stars present in 
   the narrow F502N images way saturate and are useless for our CTE 
   measurements. And further, because the master drizzled image was created 
   from F502N images, the master catalog does not contain many of F606W's 
   faint sources.  One might want to create a second master image/catalog 
   for F606W images only. This is not high priority since F606W
   observations were discontinued in 2012. But still might be nice to have
   the F606W results available for the new pipeline.
   This actually might be a little complicated when filling out the 
   photometry database. Will need to make sure that source IDs for sources
   that also appear in F502N match.

2. Speed up the database ingestion?
   I believe the slowest part to the entire pipeline is the ingestion
   of the photometry table entries.
   Since it is a SQLITE database, you cannot parallelize the ingestion
   of entries. You *could* switch to a MySQL database. I chose not to begin
   with a MySQL database because it lives on a server and I would not have
   been able to access it offline to do testing (I was working offline
   in random location A LOT while I was developing this...). It also 
   gave me the freedom to make a bunch of copies to try out one thing on
   one version of the database. In fact, currently two working versions
   exist one based on Python photometry and the other on IRAF photometry. 
   A final concern is dealing with IT to setup and maintain the MySQL 
   database (and stuff goes down a lot in this place, in spite of IT's 
   best attempts...). I like the automaty of the SQLITE database, but it's
   up to you. 
   Once IT has given you the access to MySQL, switching to MySQL is 
   literally one line of code (setting the connection to the server instead 
   of sqlite in `load_connection` function in `database_interface.py`).

3. Switch from IRAF's 'phot' routine to Python's 'aperture_photometry'
   function from the photutils package.  I already started this. The
   wrappers to the Python photometry exist in the cgosmeyer repo 
   `photutils_plus`.  I did not use the Python
   photometry in my ISRs, however because of the disreprency I found between
   the IRAF photometry (see Bajaj & Khandrika 2017), and I at the time 
   wanted consistency with the photometry from the IDL pipeline, which also 
   used the IRAF 'phot' routine. 
   In addition, I'm not 100% trusting of photutils's functions. They don't
   have all the functionality of the IRAF routines, and my wrappers to them
   in detector_tools could do only so much. For instance, centroiding from
   photutils still does not work on images with multiple sources, and I believe
   that to be a contributing factor to why my CTE measurements from Python
   photometry look so much more gross than the ones from IRAF photometry.

   I have added an option to mask cosmic rays and bad pixels on top of the 
   Python photometry (see Issues 6 and 7). CTE measurements with Python
   still look gross. 

   See Section 3.5 for how to switch between the Python and IRAF
   photometry.

4. Improve the 2-d polynomial fitting function in `uvis_external_cte_plots.py`. 
   Relevant functions are `fit_empirical_model` and `polyfit2d`.
   Right now it deals poorly with the ends, so that the low flux bin is not 
   fitting at all, and I end up tossing out that bin entirely. Probably not the 
   best way to do this. The IDL pipeline’s sfit procedure is clearly better 
   at this. (but IDL also has more sophisticated fitting procedures…).  I just 
   pulled the 2-d poly fitter off stack overflow and worked better than anything 
   else I tried.  It would be worth your time to play with this. 

5. In the 180-degree dataset, test fixing where the fit for the flux-ratio vs chip2 
   y-position crosses the 1.0 ratio to always be at 1024 pixels.
   There might be a statistical error in our measurement of CTE, which depends on 
   the fit crossing at 1024 pixels. For the long exposures this does not seem to 
   vary much image-pair to image-pair.  For the short exposures, especially the low
   flux bins, you do see the fit crossing at different y-positions, although often
   these bins are so messy anyway we don’t count them for much… I’m not 100% convinced
   this extra correction makes sense because there’s additional sources of noise
   such as cosmic rays, possibly stars being read out by different amplifiers, flat
   field fluxuations, geometric distortion, etc.  Really, might make more sense to
   quantify it as an uncertainty and report it as that, rather than try to correct
   it?  
   We do NOT do this for the nominal dataset because correcting the fit always
   cross at y-position 1024 at ratio 1.0 assumes (1) the CTE on the two chips behaves
   approximately the same and (2) the source at the center (y=1024) on one chip
   falls at the center on the second chip.  
   For the nominal monitor assumption (1) fails the CTE degredation in the two chips
   cannot be assumed to be the same because the two chips were cut from different 
   wafers and likely they had different populations of charge traps. ACS can get away
   with this assumption because the two chips on WFC were indeed cut from the same
   wafer. You can convince yourself of the inaccuracy of assumption (2) just
   by blinking the image pairs. The target fields are very often offset by many pixels.
   Note that 'fixing the fit' was not done in the original CTE monitor written by
   Kai Noeske and I believe he did this deliberately because he adopted almost all the other
   methods employed by M. Chiaberge in his external CTE monitor for ACS/WFC. 
   (see the ACS ISR 2012-05)
   But anyway you should try to see what happens when you correct the fit to always cross
   at y-position 1024 at ratio 1.0 in the 180-degree dataset. In ISR 2016-17 I gave only
   a back of the evelope calculation. 

6. Test masking bad pixels and cosmic rays. I wrote scripts to do create the masks and
   have resulting photometry and CTE measurements, but never stringently looked 
   through the cosmic ray masks to see if they weren't over-masking and never did 
   a strict comparision between the masked and non-masked photometry and CTE measurements. 
   Should also quantify how the masking affects where the linear fit crosses 1.0
   ratio (as described in Issue 5).
   The script is called `run_adriz.py` and if it is used, it needs be run as the 
   second prep step after `run_tweakreg.py`.

7. Related to Issue 6...
   If decide to make `run_adriz.py` the second prep step, might consider some of
   following changes.
   * Create seperate columns in the photometry table for bad pixel and cosmic ray-masked 
   results?
   (I already have a seperate database filled with bad pixel and cosmic-ray masked 
   Python-based photometry. This may be cleaner than creating seperate columns in the
   main database, and anyway, before investing the time, it would depend on whether 
   we think doing cosmic-ray masking is actually a good idea.)
   * Change the name of the 'MNCLIP_BKGRD' column in the FILEINFO table to 'MDRIZSKY'.
   See above in section 4.1, Procedure for Adding/Renaming a Column.

8. Make the FLT/FLC loop and flashlvl loop user-specified in `run_uvis_external_cte.py`. 
   So the user should be able to specify that they just want to run over flashlvl=6 or
   they just want to run the pipeline over FLCs. 
   It's silly I hard-coded them when everything else is set as a user-given argument
   or set in config.py.
   You'll need add them to the arg_parse function and call, to the top docstring, and 
   of course to the run_uvis_exernal_cte function.

9. Pull the hard-coded 'bad visits' from `run_uvis_external_cte.py` to some sort of
   lists in config.py, in the same vein as Issue 8.
   There are some data subsets that the pipeline is constrained from running over. 
   See the myriad of comments in `run_uvis_external_cte.py`.
   The primary subsets are the charge injection visits. From the standpoint of best
   coding practices, I should not be throwing all these data subsets out in the script
   but rather by reading some list of "forbidden visits" in config.py. It would 
   be effectively the same thing, but keeping with the idea that config.py should be the primary 
   judge of what is legal for the pipeline to run on. 
   But why not just remove the offending images entirely from the data directory?
   I thought to do this when I began setting up the pipeline, but then rejected it for 
   two primary reasons.
   (1) When copying data from MAST or Quicklook, or re-grabbing everything after re-calibrating,
   it's annoying to weed out the bad ones and humans being humans you always miss some, and
   I've done it before.  Hardcoding the forbidden lists instead of hoping I didn't copy
   charge-injected image by accident showed itself to be more dummy-proof. 
   (2) I had a vision that the pipeline should be able to run over some of the stranger
   datasets, such as the Charge Injection ones. So, being opimistic, I retained all those
   images in the data directory. The pipeline doesn't run properly on them now, but that might be
   a matter of adding some special functions to handle them (and what else is a pipeline for??). 
   The pipeline should be broad enough that it should be possible for a user to say 
   they want CTE measurements of whatever subset were observed in the 6+ years of the CTE monitor.

