#! /usr/bin/env python

"""Tests for the ``log_exoctk`` module.

Authors
-------

    Matthew Bourque

Use
---

    These tests can be run via the command line (omit the ``-s`` to
    suppress verbose output to stdout):
    ::

        pytest -s test_log_exoctk.py
"""

import os

import astropy
import sqlite3

from exoctk import log_exoctk


def test_create_db():
    """Test the ``create_db`` function"""

    # Create the database
    db_path = './test.db'
    log_exoctk.create_db(db_path)

    # Ensure the database exists
    assert os.path.exists(db_path)

    # Remove the database
    os.remove(db_path)


def test_load_db():
    """Test the ``load_db`` function"""

    # Create the database
    db_path = './test.db'
    log_exoctk.create_db(db_path)

    # Ensure the database can be loaded
    cursor = log_exoctk.load_db(db_path)
    assert isinstance(cursor, sqlite3.Cursor) and cursor is not None

    # Remove the database
    os.remove(db_path)


def test_log_form_input():
    """Test the ``log_form_input`` function"""

    # Create database
    db_path = './test.db'
    log_exoctk.create_db(db_path)

    # Define test parameters
    form_dict = {}
    tables = ['groups_integrations', 'limb_darkening', 'contam_visibility', 'phase_constraint', 'fortney', 'generic']
    database = sqlite3.connect(db_path, isolation_level=None, detect_types=sqlite3.PARSE_DECLTYPES, check_same_thread=False)
    cursor = database.cursor()

    # Test the function
    for table in tables:
        log_exoctk.log_form_input(form_dict, table, database)

        # Make sure a record was inserted
        cursor.execute(f'SELECT * FROM {table}')
        rows = cursor.fetchall()
        assert len(rows) > 0

    # Remove the database
    os.remove(db_path)


def test_scrub():
    """Test the ``scrub`` function"""

    table_name = 'DROP TABLE groups_integrations'
    returned_name = log_exoctk.scrub(table_name)

    assert returned_name == 'DROPTABLEgroupsintegrations'


def test_view_log():
    """Test the ``view_log`` function"""

    # Create database
    db_path = './test.db'
    log_exoctk.create_db(db_path)

    # Define test parameters
    tables = ['groups_integrations', 'limb_darkening', 'contam_visibility', 'phase_constraint', 'fortney', 'generic']
    conn = sqlite3.connect(db_path, isolation_level=None, detect_types=sqlite3.PARSE_DECLTYPES, check_same_thread=False)
    database = conn.cursor()

    # Test the function
    for table in tables:
        data = log_exoctk.view_log(database, table)
        assert isinstance(data, astropy.table.table.Table) and data is not None

    # Remove the database
    os.remove(db_path)
