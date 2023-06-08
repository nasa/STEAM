#! test_remove.py
#! In order to run these tests with verbose output:
#! >>> python test_remove.py -v
#!
"""
Test the capability to remove meshes, solutions, and points from a database.

Creates a database object with solutions on several meshes and then removes
some of them in various fashions.
"""

import steam
import unittest as ut
import pandas as pd
from code import interact
#interact( local = dict( globals(), **locals() ) )

class TestDatabase( ut.TestCase ):

    @classmethod
    def setUp( self ):
        """This method gets called before each of the tests is run in 
        order to create the fake data that is needed.  Here we will create 
        a database with several meshes.
        """ 

        ### Make several meshes
        m_sparse = steam.mesh.mk_cube( 1, 3 )
        m_med    = steam.mesh.mk_cube( 1, 8 )
        m_fine   = steam.mesh.mk_cube( 1, 11 )

        ### Prepare for progress bar -- 
        soln_per_grid = 5
        n_iter = 1
        total_iter = 3 * 2 * soln_per_grid

        ### Make a database and populate with 3 solutions on each mesh
        self.db = steam.database.Database()
        for m_num, mesh in enumerate([m_sparse, m_med, m_fine]):
            for store_at in ( 'ELEMENTS', 'NODES' ):
                for j in range(soln_per_grid):
                    steam.util.progress_bar( n_iter, total_iter, 
                                             'Adding solns to DB:' )
                    n_iter += 1

                    ### Make a DBPoint
                    pt = steam.database.DBPoint()
                    pt.data['orig_mesh'] = m_num        ### Original mesh
                    pt.data['vel']       = 4 + 2*j      ### Dummy velocity
                    pt.data['alpha']     = 20 + j*m_num ### Dummy alpha
                    
                    pt.mesh = mesh
                    pt.soln = steam.solution.uniform_soln( pt.mesh, store_at,
                                           values = [1.0, 'rand', 'rand'] )

                    ### Add the point to the database
                    self.db.add_point( pt )

        ### Append a few points again; some solutions shouldn't be unique
        for i, pt in enumerate( self.db ):
            self.db.add_point( pt )
            if i > 5:
                break

    def test_pt_by_ind( self ):
        """ Remove one point at a time."""
        
        ### Remove every point from the database
        for pt in self.db:
            ### Check whether the soln and mesh are unique to this point
            meshid = pt.get_meshid()
            solnid = pt.get_solnid()
            self.assertIn( meshid, self.db.meshes )
            self.assertIn( solnid, self.db.solns )
            unique_mesh = True if len( self.db.data.query(
                    'meshid == "{:s}"'.format(meshid) ) ) == 1 else False
            unique_soln = True if len( self.db.data.query(
                    'solnid == "{:s}"'.format(solnid) ) ) == 1 else False

            ### Remove the point
            self.assertIs( pt.get_database(), self.db )
            self.db.remove_pt( pt )
            self.assertNotIn( pt.get_index(), self.db.data.index )

            ### Check that the soln and mesh are missing (if they should be)
            if unique_mesh:
                self.assertNotIn( meshid, self.db.meshes )
            else:
                self.assertIn( meshid, self.db.meshes )

            if unique_soln:
                self.assertNotIn( solnid, self.db.solns )
            else:
                self.assertIn( solnid, self.db.solns )


#    def test_pts_by_ind( self ):
#        """ Remove several points at once."""
#
#        ### Remove every point in sets
#        for ind in self.db.data.index[::5]:
#            inds = list( range(ind, ind+5) )
#
#            ### Remove the points
#            self.db.remove_pt( inds )
#            self.assertTrue( 
#                all( [index not in self.db.data.index for index in inds] ) )
#
#
    def test_fake_pt( self ):
        """ Attempt to remove a point not actually in the database."""

        ### Make a point that's not in the database
        data = self.db.data.loc[0]  # Directly copy data from another point
        fake_pt = steam.database.DBPoint( data=data, 
           mesh=self.db.meshes[data.meshid], soln=self.db.solns[data.solnid] )

        self.assertRaises( AssertionError, self.db.remove_pt, fake_pt )

    @ut.skip( ' ---- Placeholder; not written yet.' )
    def test_soln( self ):
        """ Remove a solution from the database."""
        raise NotImplementedError()

    @ut.skip( ' ---- Placeholder; not written yet.' )
    def test_fake_soln( self ):
        """ Attempt to remove a solution not actually in the database."""
        raise NotImplementedError()

    @ut.skip( ' ---- Placeholder; not written yet.' )
    def test_mesh( self ):
        """ Remove a mesh from the database. """
        raise NotImplementedError()

    @ut.skip( ' ---- Placeholder; not written yet.' )
    def test_fake_mesh( self ):
        """ Attempt to remove a mesh not actually in the database."""
        raise NotImplementedError()

#! Run the tests
if __name__ == '__main__':
    ut.main()
