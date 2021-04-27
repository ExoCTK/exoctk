"""This module creates and manages a SQL database as a log for all jobs
submitted via the exoctk web app.

Authors
-------

    - Joe Filippazzo

Use
---

    This module is intended to be imported and used within a separate
    python environment, e.g.

    ::

        from exoctk import log_exoctk
        log_exoctk.create_db()

Dependencies
------------

    - ``astropy``
    - ``numpy``
    - ``pathlib``
    - ``sqlite3``
"""

import os
import datetime

import astropy.table as at
import numpy as np
from pathlib import Path
import sqlite3


def create_db(dbpath, overwrite=True):
    """Create a new database at the given ``dbpath``

    Parameters
    ----------
    dbpath : str
        The full path for the new database, including the filename
        and ``.db`` file extension.
    schema : str
        The path to the ``.sql`` schema for the database
    overwrite : bool
        Overwrite dbpath if it already exists
    """

    if dbpath != ':memory:':

        # Make sure the path is valid
        if not os.path.exists(os.path.dirname(dbpath)):
            raise IOError('Not a valid path:', dbpath)

        # Make sure the file is a .db
        if not dbpath.endswith('.db'):
            raise IOError('Please provide a path with a .db file extension')

        # Use pathlib.Path to become compliant with bandit.
        p = Path(dbpath)

        # Remove existing file if overwriting
        if os.path.isfile(dbpath) and overwrite:
            p.unlink()

        # Make the new file
        p.touch()

    # Generate the tables
    conn = sqlite3.connect(dbpath)
    cur = conn.cursor()

    # Table for groups_integrations
    cur.execute("CREATE TABLE 'groups_integrations' ('id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'targname' TEXT, 'kmag' REAL, 'mod' TEXT, 'obs_time' REAL, 'n_group' REAL, 'ins' TEXT, 'filt' TEXT, 'filt_ta' TEXT, 'subarray' TEXT, 'subarray_ta' TEXT, 'sat_mode' TEXT, 'sat_max' REAL, PRIMARY KEY(id));")

    # Table for limb_darkening
    cur.execute("CREATE TABLE 'limb_darkening' ('id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'n_bins' INTEGER, 'teff' REAL, 'logg' REAL, 'feh' REAL, 'bandpass' TEXT, 'modeldir' TEXT, 'wave_min' REAL, 'mu_min' REAL, 'wave_max' REAL, 'local_files' TEXT, 'pixels_per_bin' INTEGER, 'uniform' TEXT, 'linear' TEXT, 'quadratic' TEXT, 'squareroot' TEXT, 'logarithmic' TEXT, 'exponential' TEXT, 'three_parameter' TEXT, 'four_parameter' TEXT, PRIMARY KEY(id));")

    # Table for contam_visibility
    cur.execute("CREATE TABLE 'contam_visibility' ('id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'targname' TEXT, 'ra' REAL, 'dec' REAL, 'inst' TEXT, 'companion' TEXT, PRIMARY KEY(id));")

    # Table for phase_constraint
    cur.execute("CREATE TABLE 'phase_constraint' ('id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'targname' TEXT, 'orbital_period' REAL, 'eccentricity' REAL, 'transit_type' TEXT, 'omega' REAL, 'inclination' REAL, 'transit_time' REAL, 'window_size' REAL, 'observation_duration' REAL, 'minimum_phase' REAL, 'maximum_phase' REAL, PRIMARY KEY(id));")

    # Table for fortney grid
    cur.execute("CREATE TABLE 'fortney' ('id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'ptemp' REAL, 'pchem' TEXT, 'cloud' TEXT, 'pmass' REAL, 'm_unit' TEXT, 'refrad' REAL, 'r_unit' TEXT, 'rstar' REAL, 'rstar_unit' TEXT, PRIMARY KEY(id));")

    # Table for generic grid
    cur.execute("CREATE TABLE 'generic' ('id' INTEGER NOT NULL UNIQUE, 'date' TEXT NOT NULL, 'temperature' REAL, 'gravity' REAL, 'r_planet' REAL, 'r_star' REAL, 'condensation' TEXT, 'metallicity' REAL, 'c_o' REAL, 'haze' REAL, 'cloud' REAL, PRIMARY KEY(id));")

    # Commit and close the connection
    conn.commit()
    conn.close()

    if os.path.isfile(dbpath):
        print("ExoCTK database created at {}".format(dbpath))


def load_db(dbpath):
    """Load a database

    Parameters
    ----------
    dbpath : str
        The path to the ``.db`` database file


    Returns
    -------
    cur : ``sqlite.connection.cursor`` obj
        An SQLite3 Cursor object if dbpath is found.
    """

    if os.path.isfile(dbpath) or dbpath == ':memory:':

        con = sqlite3.connect(dbpath, isolation_level=None, detect_types=sqlite3.PARSE_DECLTYPES, check_same_thread=False)
        cur = con.cursor()

        print('Database loaded: {}'.format(dbpath))

        return cur

    else:
        print("Sorry, could not find the file '{}'".format(dbpath))


def log_form_input(form_dict, table, database):
    """A function to store the form inputs of any page ``GET`` requests
    in a database.

    Parameters
    ----------
    form_dict : dict
        The dictionary of form inputs
    table : str
        The table name to INSERT on
    database : ``sqlite.connection.cursor`` obj
        The database cursor object
    """
    try:

        # Get the column names
        colnames = np.array(database.execute("PRAGMA table_info('{}')".format(table)).fetchall()).T[1]

        # Convert hyphens to underscores and leading numerics to letters for db column names
        inpt = {k.replace('-', '_').replace('3', 'three').replace('4', 'four'): v for k, v in form_dict.items()}

        # Add a timestamp
        inpt['date'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

        # Insert the form values that are valid column names
        cols = [col for col in inpt.keys() if col in colnames]
        vals = [str(inpt.get(col, 'none')) for col in cols]
        qry = "Insert Into {} ({}) Values ({})".format(table, ', '.join(cols), ', '.join('?' * len(cols)))
        database.execute(qry, vals)

    except Exception as e:
        print('Could not log form submission.')
        print(e)


def view_log(database, table, limit=50):
    """Visually inspect the job log.

    Parameters
    ----------
    database : str or ``sqlite3.connection.cursor`` obj
        The database cursor object
    table : str
        The table name
    limit : int
        The number of records to show

    Returns
    -------
    table : ``astropy.Table`` obj
        An astropy.table object containing the results.

    """

    if isinstance(database, str):
        DB = load_db(database)
    elif isinstance(database, sqlite3.Cursor):
        DB = database
    else:
        print("Please enter the path to a .db file or a sqlite.Cursor object.")

    # Query the database
    colnames = np.array(DB.execute("PRAGMA table_info('{}')".format(table)).fetchall()).T[1]
    results = DB.execute("SELECT * FROM {} LIMIT {}".format(table, limit)).fetchall()

    # Empty table
    table = at.Table(names=colnames, dtype=['O'] * len(colnames))

    # Add the results
    if len(results) > 0:
        for row in results:
            table.add_row(row)

    return table


def scrub(table_name):
    """Snippet to prevent SQL injection attcks! PEW PEW PEW!"""
    return ''.join(chr for chr in table_name if chr.isalnum())
