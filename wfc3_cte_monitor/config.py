"""
Contains static lists, paths, and database name to be imported.
The list of the program numbers needs only be added to once a year.
"""

program_list = ['11924', '12379', '12348', '12692', '13083', '13566', '14012', '14378']

flashlvl_list = [0, 6, 12, 18, 24, 33, 55, 91, 116]

target_list = ['ngc104', 'ngc6791', 'ngc6583']

filter_list = ['F502N', 'F606W']

fluxbins_lo = [250, 500, 500, 1000, 2000, 2000, 4000, 8000]

fluxbins_hi = [500, 1000, 2000, 2000, 4000, 8000, 8000, 32000]

db_name = 'uvis_external_cte.db'

# Set path to outputs directory.
path_to_outputs = '/your/path/to/outputs/'

# Path to data directory.
path_to_data = '/your/path/to/data/'

# Clone this repo from Quicklook github and set the path here.
path_to_scripts = '/your/path/to/wfc3_cte_monitor/'

# Set path to automated_outputs in central storage.
path_to_automated_outputs = '/path/to/quicklooks/automated_outputs/'

# Path to the database.
path_to_db = '/your/path/to/database/'
