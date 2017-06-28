#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

setup(name = 'wfc3_cte_monitor',
      description = 'Hubble WFC3/UVIS CTE monitor package',
      author = 'Space Telescope Science Institute & C.M. Gosmeyer',
      url = 'https://github.com/cgosmeyer/wfc3_cte_monitor',
      packages = find_packages(),
      install_requires = ['astropy', 'numpy', 'photutils', 'pyraf',
       'scipy', 'sqlalchemy']
     )