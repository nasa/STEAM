#! /bin/env python
# test_manipulation.py
""" This set of unit tests is intended to exercise the various methods for
manipulating solutions.
"""

import unittest as ut
import steam
import numpy as np
from numpy.random import rand
import pandas as pd

cube_pts = 15

class TestManipulation( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Make a simple cube test case to work with.
        """

        ### Cubes with solution at NODES
        #   Define a mesh that DOES NOT have xyz_el populated
        cls.mesh = steam.mesh.mk_cube( 1, cube_pts)
        cls.mesh.get_xyz_el()

        ### Set up a static solution for the sphere
        cls.elem_soln = steam.solution.uniform_soln(
                           cls.mesh, 'ELEM', [1.25, 3.141] )
        cls.node_soln = steam.solution.uniform_soln(
                           cls.mesh, 'NODE', [2.01, 6.282] )

    def test_formats( self ):
        """Add a new variable to the solutions with different input formats"""

        ### New data
        node_dat = rand( self.mesh.xyz_pt.shape[0] )
        elem_dat = rand( self.mesh.xyz_el.shape[0] )

        ### Element-based solution -- add data in multiple formats
        self.elem_soln.add_var( 'list', list(elem_dat) )
        self.elem_soln.add_var( 'tuple', tuple(elem_dat) )
        self.elem_soln.add_var( 'array', np.array(elem_dat) )
        self.elem_soln.add_var( 'series', pd.Series(elem_dat) )

        ### Check that all added columns are equal
        self.assertTrue( self.elem_soln.data.list.equals( 
                         self.elem_soln.data.tuple) )
        self.assertTrue( self.elem_soln.data.tuple.equals( 
                         self.elem_soln.data.array) )
        self.assertTrue( self.elem_soln.data.array.equals( 
                         self.elem_soln.data.series) )
        self.assertTrue( self.elem_soln.data.series.equals( 
                         self.elem_soln.data.list) )

        ### Do the same for Node-based solution
        self.node_soln.add_var( 'list', list(node_dat) )
        self.node_soln.add_var( 'tuple', tuple(node_dat) )
        self.node_soln.add_var( 'array', np.array(node_dat) )
        self.node_soln.add_var( 'series', pd.Series(node_dat) )

        ### Check that all added columns are equal
        self.assertTrue( self.node_soln.data.list.equals( 
                         self.node_soln.data.tuple) )
        self.assertTrue( self.node_soln.data.tuple.equals( 
                         self.node_soln.data.array) )
        self.assertTrue( self.node_soln.data.array.equals( 
                         self.node_soln.data.series) )
        self.assertTrue( self.node_soln.data.series.equals( 
                         self.node_soln.data.list) )

    def test_add_with_index( self ):
        """ Use add_var with option index_in input """

        ### New data
        elem_dat = rand( self.mesh.xyz_el.shape[0] )
        
        ### Add data with the same index as the solution
        self.elem_soln.add_var( 'full_dat', elem_dat, 
                                index_in=self.elem_soln.data.index )
        for val in self.elem_soln.data.full_dat:
            self.assertFalse( np.isnan( val ) )
        
        ### Add discontinuous data and make sure it appears in the right place
        disc_dat = list( range(0,self.node_soln.data.shape[0],2) )
        self.node_soln.add_var( 'discontinuous', disc_dat, index_in=disc_dat )
        for ind, val in self.node_soln.data.discontinuous.items():
            if ind in disc_dat:
                self.assertEqual( val, ind )
            else:
                self.assertTrue( np.isnan( val ) )

    def test_redundant( self ):
        """ Add the same column twice and ensure an error is raised """

        ### Create and add the data
        node_dat = rand( self.mesh.xyz_pt.shape[0] )
        self.node_soln.add_var( 'data', node_dat )

        with self.assertRaises( ValueError ) as cm:
            self.node_soln.add_var( 'data', node_dat+1.0 )

#! Run the tests
if __name__ == '__main__':
    ut.main()

