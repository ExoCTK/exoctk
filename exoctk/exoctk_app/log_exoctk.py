"""
This module creates and manages a SQL database as a log for all jobs submitted via the exoctk web app
"""
import os
import numpy as np
import sqlite3
import datetime
import astropy.table as at

def create_db(dbpath, overwrite=True):
    """
    Create a new database at the given dbpath

    Parameters
    ----------
    dbpath: str
        The full path for the new database, including the filename and .db file extension.
    schema: str
        The path to the .sql schema for the database
    overwrite: bool
        Overwrite dbpath if it already exists
    """
    # Make sure the path is valid
    if not os.path.exists(os.path.dirname(dbpath)):
        raise IOError('Not a valid path:',dbpath)
        
    # Make sure the file is a .db
    if not dbpath.endswith('.db'):
        raise IOError('Please provide a path and file name with a .db file extension')
        
    # Remove existing file if overwriting
    if os.path.isfile(dbpath) and overwrite:
        os.system('rm {}'.format(dbpath))
    
    # Make the new file
    os.system('touch {}'.format(dbpath))
    
    # Generate the tables
    conn = sqlite3.connect(dbpath)
    cur = conn.cursor()
    cur.execute("CREATE TABLE 'tor' ( 'id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'ins' TEXT, 'mag' REAL, 'groups' TEXT, 'amps' INTEGER, 'subarray' TEXT, 'sat_lvl' REAL, 'sat' TEXT, 'T' REAL, 'n_reset' INTEGER, 'band' TEXT, 'filt' TEXT, PRIMARY KEY(id));")
    cur.execute("CREATE TABLE 'ldc' ( 'id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'n_bins' INTEGER, 'teff' REAL, 'logg' REAL, 'feh' REAL, 'bandpass' TEXT, 'modeldir' TEXT, 'wave_min' REAL, 'mu_min' REAL, 'wave_max' REAL, 'local_files' TEXT, 'pixels_per_bin' INTEGER, 'uniform' TEXT, 'linear' TEXT, 'quadratic' TEXT, 'squareroot' TEXT, 'logarithmic' TEXT, 'exponential' TEXT, 'three_parameter' TEXT, 'four_parameter' TEXT, PRIMARY KEY(id));")
    conn.close()
    
    if os.path.isfile(dbpath):
        print("ExoCTK database created at {}".format(dbpath))

def load_db(dbpath):
    """
    Load a database
    
    Parameters
    ----------
    dbpath: str
        The path to the .db database file
    """
    if os.path.isfile(dbpath):
    
        con = sqlite3.connect(dbpath, isolation_level=None, detect_types=sqlite3.PARSE_DECLTYPES, check_same_thread=False)
        cur = con.cursor()
        
        return cur
        
    else:
        print("Sorry, could not find the file '{}'".format(dbpath))

def log_form_input(form_dict, table, database):
    """
    A function to store the form inputs of any page GET requests in a database
    
    Parameters
    ----------
    form_dict: dict
        The dictionary of form inputs
    table: str
        The table name to INSERT on
    database: sqlite.connection.cursor
        The database cursor object
    """
    try:
        # Convert hyphens to underscores and leading numerics to letters for db column names
        inpt = {k.replace('-','_').replace('3','three').replace('4','four'):v[0] for k,v in dict(form_dict).items()}
    except TypeError:
        inpt = form_dict
    
    # Add a timestamp
    inpt['date'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    
    # Insert the form valsues
    qmarks = ', '.join('?' * len(inpt))
    qry = "Insert Into {} ({}) Values ({})".format(table, ', '.join(inpt.keys()), qmarks)
    
    database.execute(qry, list(inpt.values()))
    
def view_log(database, table):
    """
    Visually inspect the job log
    
    Parameters
    ----------
    database: str, sqlite3.connection.cursor
        The database cursor object
    table: str
        The table name
    """
    if isinstance(database, str):
        DB = load_db(database)
    elif isinstance(database, sqlite3.Cursor):
        DB = database
    else:
        print("Please enter the path to a .db file or a sqlite.Cursor object.")
        
    # Query the database
    table = scrub(table)
    colnames = np.array(database.execute("PRAGMA table_info('{}')".format(table)).fetchall()).T[1]
    results = database.execute("SELECT * FROM {}".format(table)).fetchall()
    
    results = at.Table(list(np.array(results).T), names=colnames)
    
    # Print them
    return results
    
def scrub(table_name):
    """
    Snippet to prevent SQL injection attcks! PEW PEW PEW!
    """
    return ''.join( chr for chr in table_name if chr.isalnum() )
    