#!/usr/bin/env python
### test_table.py
###
### Test the basic functionality of the Table class.
### In order to run these tests with verbose output:
### >>> python test_table.py -v
###
### NOTE: The assertCountEqual only works in unittest version 3.2 and 
###       therefore will not work in the python 2.7 environment.  Therefore,
###       tests that rely on it will be skipped.

import scipy as sp
import unittest as ut

from steam import table
from code import interact

class TestTable( ut.TestCase ):

    ##########################################
    def setUp( self ):
        """This method gets called before each test_* method in order to prepare
        the proper variables/data.  Therefore, it will initialize the two 
        tables that will be added.
        """ 

        ### Initialize and populate the tables
        self.tab1 = table.Table()
        self.tab1.read_cdat( 'table1.dat' )

        self.tab2 = table.Table()
        self.tab2.read_cdat( 'table2.dat' )

        self.tab3 = table.Table()
        self.tab3.read_cdat( 'table1.dat' )  # initially identical to self.tab1

    ##########################################
    @ut.skipUnless( 'assertCountEqual' in dir( ut.TestCase ), 
                    'Running an old version of '
                    'unittest that does not have TestCase.assertCountEqual' )
    def test_read_cdat( self ):
        """Call the read_cdat method and check the results. """
        
        self.assertCountEqual( self.tab1.data.columns, self.tab2.data.columns )

    ##########################################
    def test_compare_varlists( self ):
        """Check the ability to compare the variable lists of tables. """

        self.assertTrue( self.tab1.compare_varlists( self.tab1 ) )
        self.assertTrue( self.tab1.compare_varlists( self.tab2 ) )

    ##########################################
    def test_add_table( self ):
        """Add self.tab2 to self.tab3 and check the results. """
    
        self.tab3.add_table( self.tab2 )

        self.assertEqual( self.tab3.data.shape[1], self.tab1.data.shape[1] )
        self.assertEqual( self.tab3.data.shape[0], 
                          self.tab1.data.shape[0] + self.tab2.data.shape[0] )

    ##########################################
    def test_bad_add( self ):
        """Add incompatible tables and check the results. """
        
        ### Set up tab3 with bad data
        self.tab3 = table.Table()
        self.tab3.read_cdat( 'table2_bad.dat' )

        self.assertRaises( ValueError, self.tab1.add_table, self.tab3 )

    ##########################################
    def test_add_vars( self ):
        """Add an arbitrary variable column to tab1. """

        ### Add one variable
        self.tab1.add_vars( ['Test1'], sp.arange( self.tab1.data.shape[0]
                                       ).reshape( self.tab1.data.shape[0],1) )
        self.assertEqual( self.tab1.data.shape[1], self.tab2.data.shape[1] + 1)
        self.assertIn( 'Test1', self.tab1.data.columns )

        ### Add 2 variables simultaneously
        self.tab1.add_vars( ['Test2', 'Test3'], self.tab1.data.values[:,:2] )
        self.assertEqual( self.tab1.data.shape[1], self.tab2.data.shape[1] + 3)
        self.assertIn( 'Test2', self.tab1.data.columns )
        self.assertIn( 'Test3', self.tab1.data.columns )

    ##########################################
    def test_remove_vars( self ):
        """Remove a real variable from tab1 and attempt to remove fake vars.

        """

        ### Remove variables
        remove_vars = ['Alt(km)']
        deleted = self.tab1.remove_vars( var_names = remove_vars )
        for var in remove_vars:
            self.assertNotIn( var, self.tab1.data.columns )

        ### Attempt to remove variables that aren't present
        missing_vars = ['missing', 'variables']
        self.assertRaises(KeyError, self.tab1.remove_vars, missing_vars)

### Run the tests
if __name__ == '__main__':
    ut.main()
