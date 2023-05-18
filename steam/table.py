### Table.py
###
### Cameron Ebner
### NASA-JSC-EG3
### cameron.j.ebner@nasa.gov
###
### This file contains the class definition for the "Table" object within
### EG3's as-yet unnamed database tool.
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import numpy as np
import scipy as sp
import pandas as pd
import steam.util

#from code import interact
#interact( local = dict( globals(), **locals() ) )

################################################################################
class Table:
    """The table object is intended to wrap a pandas DataFrame object and
    provide functionality that will be useful for users of the STEAM tool.
    These capabilities will include the abilities to read CSV tables from text
    files, calculate new table variables, add multiple tables together, and
    store the tables in either text or HDF5 format.
    """

    ##################################################
    def __init__( self, data=None):

        """Initialize the Table attributes.

        The CDAT file to read is the meat of the Table.  It is a 2-D array containing all the independent variable data for each case.

        Args:
          data (:obj:`string` or :obj:`~pandas.DataFrame`): CDAT File to read or Dataframe to inherit
        """

        ### Initialize variables
        self.data      = None           # DataFrame - Indep/Dep variables

        self.indep     = []             # Independent variables
        self.dep       = []             #   Dependent variables

        self.hdf5_file = None           # File name of HDF5 file
        self.hdf5_root = ""             # Database root in HDF5 file
        self.hdf5      = None           # File pointer for HDF5 file

        self.meta      = ""             # Metadata or description

        from os.path import isfile
        if data is None:
            pass
        elif isinstance(data, pd.DataFrame):
            self.data = data
        elif isfile(data):
            self.read_cdat(data)
        else:
            raise TypeError( '"data" input should be either a pandas.Dataframe'
                             + ' or a path to a CDAT file for reading.' )

    ##################################################
    def __str__( self ):
        """Return the print string for a Table instance.
        """
        return self.data.__str__()

    def __eq__( self, other ):
        """ Used for equality determination.
        """

        return (
                self.__class__ == other.__class__ and
                self.data.equals( other.data )
               )

    def __len__( self ):
        """ Return the number of rows in self.data"""
        return len(self.data)

    def copy(self,copy_from):
        """ Copy the data in 'copy_from' to self.
        
        NOTE: As of 12/11/18, nothing is written in this method.  I, Cam, just
              added the NotImplementedError and moved on with my life.
        """

        from copy import deepcopy

        #! This is silly.  But I don't know how else to copy the values

        raise NotImplementedError("This method doesn't actually exist. Sorry.")

    ##################################################
    def compare_varlists( self, new_table ):
        """Given another Table instance, compare the column names
        from each of them.  If the columns match perfectly, return True.  
        Otherwise, return False.

        Args:
            new_table (:obj:`~pandas.Series` or :obj:`~pandas.DataFrame`): Item to compare to current data table.

        Returns:
            :obj:`bool`: True if columns are identical between current table and new_table.
                False otherwise.
        """

        new_keys = None
        if ( isinstance(new_table,pd.Series) or 
            isinstance(new_table,pd.DataFrame) ):
            new_keys = new_table.keys()
        elif isinstance(new_table,Table):
            new_keys = new_table.data.columns
        else:
            raise Exception(
                "Attemping to add unsupported table type " +
                "to existing steam.Table!\n"
                )

        #print( '\ntype(new_keys):', type( new_keys), '\n' )
        return len( self.data.columns.difference( new_keys ) ) == 0

    ##################################################
    def add_table( self, new_table ):
        
        """Given a Table, append it to the current Table if
        the variable lists are identical.

        Args:
            new_table (:obj:`steam.table`): Table to append to current table.
        """

        ### Verify that each column in the two Tables matches
        if self.compare_varlists( new_table ):
            self.data = self.data.append( new_table.data, ignore_index=True,
                                          sort=True )
        else:
            raise ValueError( 'These classes are incompatible; Check out ' + 
                              'their columns:\n' + 
                              str( self.data.columns ) + '\n' + 
                              str( new_table.data.columns ) )

    ##################################################
    def add_vars( self, new_names, new_values ):

        """Add one or more variables to the already-populated Table.
        This is done by creating a new DataFrame from the new variables and
        using the combine_first method with the current data.

        The new data must be provided in a 2-D array with an equal number of 
        rows as data.  New variable names should be a 1-D list or array
        with one entry per column of new data.

        Args:
            new_names (:obj:`list`): List of headers for new columns.
            new_values (:obj:`numpy.array`): 2-D array of data to populate columns.
        """
        
        if len( new_names ) != new_values.shape[1]:
            raise Exception( 'New variable name(s) must be equal in length ' +
                             'to the number of columns of new data' )
        elif sp.array([i in self.data.columns for i in new_names]).any():
            raise Exception( 
                        'New variable name(s) already in table.data.columns' )

        new_frame = pd.DataFrame( new_values, columns = new_names, 
            index = self.data.index )

        ### Note: overwriting is prevented by using the combine_first method,
        ###       which prioritizes the old data over new
        self.data = self.data.combine_first( new_frame )

    ##################################################
    def remove_vars( self, var_names ):

        """Remove one or more variables from the populated Table.
        This can is done by supplying variable name(s); indices are not
        very meaningful when using DataFrames.

        Return the names of removed variables.  They returned list should
        be identical to the list passed in.

        Args:
            var_names (:obj:`list`): List of variables to remove from data table.
        Returns:
            :obj:`list`: List of removed variables
        """
        
        try:
            self.data = self.data.drop( var_names, axis = 1 )
        except:
            ### Catch any variable names that aren't actually in the data
            bad_vars = [i for i in var_names if i not in self.data.columns ]
            print( 'Attempting to delete variable(s) not ' + 
                   'stored in data!\nBad Variables: ' + str( bad_vars ) + 
                   '\nData Variables: ' +  str( self.data.columns ) )
            raise

        return

    ##################################################
    ### The following methods were stolen and adapted from Alan's 
    ### ATDatabase class
    ##################################################

    ### ## Data input
    def read_cdat(self,file):
        """Read a columnized datafile.

        This assumes that the file is space delimited and has
        one header row.  Example::
        
            Head1 Head2 Head3 ...
            val11 val12 val13 ...
            ...   ...   ...

        The header row may be commented, but has to be immediately
        preceeding the data.  Another example::

        
            # Head1 Head2 Head3 ...
            val11 val12 val13 ...
            ...   ...   ...


        Args:
            file (:obj:`string`): Path to file to be read.
        """

        # Assume that the header row is commented and try to read headers
        headers = []
        for line in open(file):
            if line[:1] == "#":
                headers = line.replace('#','').split()
                continue

            # These are only the headers if the next line does NOT have a comment
            # And has the same number of elements
            if len(headers) > 0:
                sline = line.split()
                if len(sline) == len(headers):
                    # Done!
                    break
                else:
                    headers = []
                    break
        # print(headers)

        if len(headers)> 0:
            self.data = pd.read_csv(file,sep='\s+',comment='#', names=headers,header=None)
            return

        # If I didn't find headers, then assume the header line was not commented.

        ### TODO : Add check on filetype and add as input
        self.data = pd.read_table(file,delim_whitespace=True,comment='#')

        ### Error check on number of columns
        if (len(self.data.columns) <= 0):
            raise IOError(
                " * cdat file '{}' appears to have no columns".format( file ))

    ### ## Variable Management Methods

    ##################################################
    def set_indep(self,vars):
        """Set the independent variables for the database.

        This needs to be a subset of the varaibles that are already
        defined.  If not, then it will throw a warning.  This is
        used for creating the interpolation stencils later.
        If independent variables are already set, overwrite them.
        
        Args:
            vars (:obj:`list`): Verbose list of all independent variables
        """

        ### Check to make sure that the vars are in the dataframe
        if (not steam.util.vars_in_list(vars,list(self.data.columns))):
            print("Vars: ",vars)
            print("Database: ",self.data.columns)
            raise Exception("Specified dep variables not in database!")
            return

        self.indep = list(vars)

        return

    ##################################################
    def set_dep(self,vars):
        """Set the dependent variables for the database.

        These are the variables that will be returned during
        interpolation.

        This needs to be a subset of the varaibles that are already
        defined.  If not, then it will throw a warning.  This is
        used for creating the interpolation stencils later.
        If independent variables are already set, overwrite them.
        
        Args:
            vars (:obj:`list`): Verbose list of all dependent variables
        """

        ### Check to make sure that the vars are in the dataframe
        if (not steam.util.vars_in_list(vars,list(self.data.columns))):
            print("Vars: ",vars)
            print("Database: ",self.data.columns)
            raise Exception("Specified indep variables not in database!")
            return

        self.dep = list(vars)

        return

    ### ## HDF5 Methods

    ##################################################
    def write_hdf5(self,hdf5,root="/"):
        """ Simple routine to write all table to HDF5 file.
        
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`string`,optional): Location in HDF5 to write table.
                Default is to write to root location "/".
        """

        ### Write this table to disk.  Currently, it is made
        ### queryable.
        hdf5.put(
                root+"/data",
                self.data,
                format='t',
                )
        return

    ##################################################
    def read_hdf5(self,hdf5,root="/"):
        """ Read the data in HDF5 file to memory.
        
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`string`,optional): Location in HDF5 to write table.
                Default is to read from root location "/".
        """

        self.data = hdf5.get(root+"/data")
        return
