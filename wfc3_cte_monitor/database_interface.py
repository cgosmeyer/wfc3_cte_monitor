"""
Module to interface and connect to the WFC3/UVIS external CTE monitor's 
database.  When run as script on command line, creates initial empty
database.

Each target cluster in this monitor has its own four tables in the 
database containing (1) file information, (2) photometry outputs, (3) CTE 
slope results, and (4) master catalog of known sources'  RAs and Decs 
from a 'truth' drizzled image.

M. Bourque explains the magic best:

This module serves as the interface and connection module to the
database.  The ``load_connection()`` function within allows the user
to conenct to the database via the ``session``, ``base``, and
``engine`` objects (described below).  The classes within serve as the
object-relational mappings (ORMs) that define the individual tables of
the database, and are used to build the tables via the ``base`` object.

The ``engine`` object serves as the low-level database API and perhaps
most importantly contains dialects which allows the sqlalchemy module
to communicate with the database.

The ``base`` object serves as a base class for class definitions.  It
produces ``Table`` objects and constructs ORMs.

The ``session`` object manages operations on ORM-mapped objects, as
construced by the ``base``.  These operations include querying, for
example.


Authors:

    C.M. Gosmeyer (2016)
    M. Bourque

Use:

    Run this as a script in order to create an initial database with all
    empty tables (see the Main).

    ::

    >>> python database_interface.py    

    But primarily this module is intended to be imported from various 
    uvis_external_cte modules and scripts in order for the user to interact
    with the already-built database. These are the importables objects:

    ::

    from database_interface import get_session
    from database_interface import load_connection

    from database_interface import NGC104_FileInfo
    from database_interface import NGC104_Master
    from database_interface import NGC104_Phot
    from database_interface import NGC104_Results

    from database_interface import NGC6583_FileInfo
    from database_interface import NGC6583_Master
    from database_interface import NGC6583_Phot
    from database_interface import NGC6583_Results

    from database_interface import NGC6791_FileInfo
    from database_interface import NGC6791_Master
    from database_interface import NGC6791_Phot
    from database_interface import NGC6791_Results

References:

    http://zetcode.com/db/sqlalchemy/orm/
    http://pythoncentral.io/introductory-tutorial-python-sqlalchemy/

Notes:

    Functions and documentation heavily emulated from M. Bourque's
    hstlc ``database_interface.py`` and various ``pyql`` functions.

"""

import os

from sqlalchemy import create_engine
from sqlalchemy import Boolean
from sqlalchemy import Column
from sqlalchemy import Date
from sqlalchemy import DateTime
from sqlalchemy import ForeignKey
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy import Float
from sqlalchemy import Date
from sqlalchemy import Time
from sqlalchemy import UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import relationship

from config import db_name
from config import path_to_db


__all__ = ['NGC104_FileInfo', 'NGC104_Master', 'NGC104_Phot', 'NGC104_Results', 
           'NGC6791_FileInfo', 'NGC6791_Master', 'NGC6791_Phot', 'NGC6791_Results',
           'NGC6583_FileInfo', 'NGC6583_Master', 'NGC6583_Phot', 'NGC6583_Results']

#-------------------------------------------------------------------------------# 

def get_session():
    """Return the ``session`` object of the database connection

    In many cases, all that is needed is the ``session`` object to
    interact with the database.  This function can be used just to
    establish a connection and retreive the ``session`` object.

    Returns:
        session : sqlalchemy.orm.session.Session
            Provides a holding zone for all objects loaded or associated
            with the database.

    Notes:
        From M. Bourque's 'database_interface.py'
    """

    session, base, engine = load_connection(db_name)

    return session

#-------------------------------------------------------------------------------# 

def load_connection(db_name, echo=False):
    """Create and return a connection to the database given in the
    connection string.

    Parameters:
        db_name : string
            Name of the SQL Alchemy database in '<name>.db' format.
        echo : bool
            Show all SQL produced

    Returns:
        session : sesson object
            Provides a holding zone for all objects loaded or 
            associated with the database.
        base : base object
            Provides a base class for declarative class definitions.
        engine : engine object
            Provides a source of database connectivity and behavior.

    Notes:
        From M. Bourque's 'database_interface.py'

    """
    #Four slashes for absolute path! Three slashes for relative path.
    #engine = create_engine('sqlite:///{}'.format(db_name) , echo=echo)
    engine = create_engine('sqlite:////{}/{}'.format(
        path_to_db, db_name), echo=echo)
    base = declarative_base(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    return session, base, engine

session, base, engine = load_connection(db_name)


#-------------------------------------------------------------------------------# 
# Define ORMS for Master, Phot, and File Info tables.
#-------------------------------------------------------------------------------# 

def define_columns(class_attributes_dict, header_name):
    """Dynamically define the class attributes for the ORMs.

    Parameters:
        class_attributes_dict : dictionary
            Don't worry about it.
        header_name : string
            Rootname of the column definitions text file, located in 
            'table_definitions' subdirectory.

    Returns:
        class_attributes_dict : dictionary
            Don't worry about it.
            
    Outputs:
        nothing

    Notes:
        From pyql.database.databse_interface (fall 2015 commits).
    """
    with open(os.path.join(os.path.split(__file__)[0], 'table_definitions',
                            header_name + '.txt'), 'r') as f:
        data = f.readlines()
    header_keyword_list = [item.strip().split(', ') for item in data]
    for header_keyword in header_keyword_list:
        if header_keyword[1] == 'Integer':
            class_attributes_dict[header_keyword[0].lower()] = Column(Integer())
        elif header_keyword[1] == 'String':
            class_attributes_dict[header_keyword[0].lower()] = Column(String(100))
        elif header_keyword[1] == 'Float':
            class_attributes_dict[header_keyword[0].lower()] = Column(Float())
        elif header_keyword[1] == 'Date':
            class_attributes_dict[header_keyword[0].lower()] = Column(Date())
        elif header_keyword[1] == 'Time':
            class_attributes_dict[header_keyword[0].lower()] = Column(Time())
        elif header_keyword[1] == 'DateTime':
            class_attributes_dict[header_keyword[0].lower()] = Column(DateTime)
        else:
            raise ValueError('header keyword type not recognized: {}:{}'.format(
                header_keyword[0], header_keyword[1]))

    return class_attributes_dict


#-------------------------------------------------------------------------------# 

def orm_fileinfo_factory(class_name, master_table, master_name):
    """
    Creates SQLA ORM Classes.

    Parameters:
        class_name : string
            Name of the FileInfo class table.
        master_table : Master table class
            Corresponding Master table class, so can link 'id' column.
        master_name : string
            Name of the corresponding Master table.

    Returns:
        tuple of attributes    
    """

    class_attributes_dict = {}
    class_attributes_dict['id'] = Column(Integer, primary_key=True, index=True)
    class_attributes_dict['__tablename__'] = class_name.lower()
    class_attributes_dict = define_columns(class_attributes_dict, class_name.split('_')[1])

    return type(class_name.upper(), (base,), class_attributes_dict)


#-------------------------------------------------------------------------------# 

def orm_master_factory(class_name):
    """ Creates master classes.

    Parameters:
        class_name : string
            Name of the Master table class.

    Returns:
        tuple of attributes
    """
    class_attributes_dict = {}
    class_attributes_dict['id'] = Column(Integer, primary_key=True, index=True)
    class_attributes_dict['__tablename__'] = class_name.lower()
    class_attributes_dict = define_columns(class_attributes_dict, class_name.split('_')[1])

    return type(class_name.upper(), (base,), class_attributes_dict)


#-------------------------------------------------------------------------------# 

def orm_phot_factory(class_name, master_table, master_name, 
                        fileinfo_table, fileinfo_name):
    """
    Creates SQLA ORM Classes.

    Parameters:
        class_name : string
            Name of the FileInfo class table.
        master_table : Master table class
            Corresponding Master table class, so can link 'id' column.
        master_name : string
            Name of the corresponding Master table.
        fileinfo_table : FileInfo table class
            Corresponding FileInfo table class, so can link 'imagename' column.
        fileinfo_name : string
            Name of the corresponding FileInfo table.

    Returns:
        tuple of attributes   
    """

    class_attributes_dict = {}
    class_attributes_dict['id'] = Column(Integer, primary_key=True, index=True)
    class_attributes_dict['master_id'] = Column(
        Integer, ForeignKey('{}.id'.format(master_name)), nullable=False, index=True)
    class_attributes_dict['imagename'] = Column(
        Integer, ForeignKey('{}.imagename'.format(fileinfo_name)), nullable=False, index=True)
    class_attributes_dict['__tablename__'] = class_name.lower()
    class_attributes_dict = define_columns(class_attributes_dict, class_name.split('_')[1])

    return type(class_name.upper(), (base,), class_attributes_dict)


#-------------------------------------------------------------------------------# 

def orm_results_factory(class_name):
    """
    Creates SQLA ORM Classes.

     Parameters:
        class_name : string
            Name of the Results class table.

    Returns:
        tuple of attributes      
    """

    class_attributes_dict = {}
    class_attributes_dict['id'] = Column(Integer, primary_key=True, index=True)
    #class_attributes_dict['imagename'] = Column(
    #    Integer, ForeignKey('{}.imagename'.format(fileinfo_name)), nullable=False, index=True)
    class_attributes_dict['__tablename__'] = class_name.lower()
    class_attributes_dict = define_columns(class_attributes_dict, class_name.split('_')[1])

    return type(class_name.upper(), (base,), class_attributes_dict)


#-------------------------------------------------------------------------------# 
# Initialize classes.
#-------------------------------------------------------------------------------# 

NGC104_Master = orm_master_factory('ngc104_master')
NGC6583_Master = orm_master_factory('ngc6583_master')
NGC6791_Master = orm_master_factory('ngc6791_master')

NGC104_FileInfo = orm_fileinfo_factory('ngc104_fileinfo', NGC104_Master, 
                                       'ngc104_master')
NGC6583_FileInfo = orm_fileinfo_factory('ngc6583_fileinfo', NGC6583_Master, 
                                         'ngc6583_master')
NGC6791_FileInfo = orm_fileinfo_factory('ngc6791_fileinfo', NGC6791_Master, 
                                        'ngc6791_master')

NGC104_Phot = orm_phot_factory('ngc104_phot', NGC104_Master, 
                               'ngc104_master', NGC104_FileInfo,
                               'ngc104_fileinfo')
NGC6583_Phot = orm_phot_factory('ngc6583_phot', NGC6583_Master, 
                                'ngc6583_master', NGC6583_FileInfo,
                                'ngc6583_fileinfo')
NGC6791_Phot = orm_phot_factory('ngc6791_phot', NGC6791_Master, 
                                'ngc6791_master', NGC6791_FileInfo,
                                'ngc6791_fileinfo')

NGC104_Results = orm_results_factory('ngc104_results')
NGC6583_Results = orm_results_factory('ngc6583_results')
NGC6791_Results = orm_results_factory('ngc6791_results')

#-------------------------------------------------------------------------------# 
# Main to create initial empty database
#-------------------------------------------------------------------------------# 

def create_database(base):
    """
    Creates initial database with all empty tables.
    To check that it was created properly, 

    % sqlite3 uvis_external_cte.db
    sqlite> .tables
    (Should see all tables.)
    sqlite> pragma table_info(ngc104_fileinfo);
    (Should see all columns. Check each table that they all contain wanted columns)

    Parameters:
        base : base object
            Don't worry about it.

    Returns:
        nothing

    Outputs:
        Initial 'uvis_external_cte.db' with empty tables.
    """
    base.metadata.create_all()


#-------------------------------------------------------------------------------# 

if __name__=="__main__":

    create_database(base)

