#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
A module for managing references in ExoCTK
"""
import os
import pkg_resources

import bibtexparser as bt


class References(object):
    """
    Creates and manages a References object to track references
    within an ExoCTK user session

    Attributes
    ----------
    bibfile: str
        The path to the bibtex file from which the references will be read
    refs: list
        The list of bibcodes saved during the user session
    database: bibtexparser.bibdatabase.BibDatabase object
        The database of parsed bibtex entries
    bibcodes: list
        The list of all bibcodes in the database

    """
    def __init__(self, bibfile=''):
        """
        Initializes an empty References object which points to a
        .bib file

        Parameters
        ----------
        bibfile: str
          The path to the bibtex file from which the references will be read

        """
        bibfile = bibfile or \
            pkg_resources.resource_filename('ExoCTK', 'data/core/bibtex.bib')

        # Attributes for the filepath and references
        self.bibfile = bibfile
        self.refs = []

        # Load the bibtex into a database
        bf = open(bibfile)
        self.database = bt.load(bf)
        bf.close()

        # The list of all bibcodes in the bibfile
        self.bibcodes = [i['ID'] for i in self.database.entries]

    def add(self, bibcode):
        """
        Adds a bibcode to the References object

        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be added

        """
        # Check that the bibcode is in the bibtex file
        if bibcode in self.bibcodes:
            self.refs += [bibcode]
            print(bibcode, 'added to list of references.')

        # Suggest adding it to the bibfile
        else:
            print(bibcode, 'not in bibfile at', self.bibfile)
            print('Add the bibtex entry to the file and try agin.')

    def remove(self, bibcode):
        """
        Removes a bibcode from the References object

        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be removed

        """
        # Check that the bibcode is in the bibtex file
        if bibcode in self.bibcodes:
            self.refs = [r for r in self.refs if r != bibcode]
            print(bibcode, 'removed from list of references.')

        # Nothing to remove!
        else:
            print(bibcode, 'not in bibfile at', self.bibfile)

    def write(self, filepath=''):
        """
        Write the .bib file

        Parameters
        ----------
        filepath: str
            The directory to which the .bib file should be written.
            Can be an existing .bib file or just a directory

        """
        # Use existing .bib file or create new one
        if filepath.endswith('.bib'):
            bibfile = filepath
        else:
            bibfile = os.path.join(filepath, 'biblio.bib')

        # Create a new database instance
        final = bt.bibdatabase.BibDatabase()

        # Get the relevant bibtex entries
        final.entries = [d for d in self.database.entries
                         if d['ID'] in list(set(self.refs))]

        # Write the bibtex to file
        with open(bibfile, 'w') as out:
            out.write(bt.bwriter.BibTexWriter().write(final))
