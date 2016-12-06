class References(object):
    def __init__(self, bibfile=''):
        """
        Creates a References object to track references within an ExoCTK session
    
        Parameters
        ----------
        bibfile: str
          The path to the bibtex file from which the references will be read
    
        """
        # Attributes for the filepath, file contents, and references
        self.bibfile = bibfile
        self.bibtex = pickle.load(open(bibfile, 'rb'))
        self.refs = []

    def add(self, bibcode):
        """
        Adds a bibcode to the References object
        
        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be added
        
        """
        # Check that the bibcode is i the bibtex file
        if bibcode in self.bibtex:
            self.refs += bibcode
            print(bibcode,'added to list of references.')
        
        # Suggest adding it to the bibfile
        else:
            print(bibcode,'not in bibfile at',self.bibfile)
            print('Add the bibtex entry to',self.bibfile,'and try agin.')
            
    def remove(self, bibcode):
        """
        Removes a bibcode from the References object
        
        Parameters
        ----------
        bibcode: str
            The unique compact identifier for the reference to be removed
        
        """
        # Check that the bibcode is i the bibtex file
        if bibcode in self.bibtex:
            self.refs = [r for r in self.refs if r!=bibcode]
            print(bibcode,'removed from list of references.')
        
        # Suggest adding it to the bibfile
        else:
            print(bibcode,'not in bibfile at',self.bibfile)
                
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
        bibfile = filepath if filepath.endswith('.bib') else filepath+'biblio.bib'
        
        # Iterate through the references and write the relevant bibtex to file
        for bibcode in self.refs:
            with open(bibfile, 'a') as bib:
                bib.write(self.bibtex[bibcode])
