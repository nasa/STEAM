#! test_to_common_mesh.py
#! In order to run these tests with verbose output:
#! >>> python test_to_common_mesh.py -v
#!
"""
Test the capability to convert solutions on multiple meshes in a database 
to a common mesh.

Creates a database object with solutions on several meshes and then converts
them to a different mesh.
"""

import steam
import unittest as ut
import pandas as pd
#from code import interact

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
        soln_per_grid = 10
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

        ### Check for a LibMeshInit object and delete it if necessary
        try:
            steam.libmesh._lmp
        except AttributeError:
            ### There is no Init object, which is what we want
            pass
        else:
            del(steam.libmesh._lmp)
 
    def test_raises_runtime( self ):
        """ to_common_mesh should raise a RuntimeError if LibMesh is loaded"""
        m_target   = steam.mesh.mk_cube( 1, 11 )

        ### When a LibMeshInit exists it's stored in steam.libmesh._lmp
        try:
            steam.libmesh._lmp
        except AttributeError:
            if steam.has_libmesh:
                ### Use a real LibMeshInit
                steam.libmesh.Init()
            else:
                ### Fake it
                steam.libmesh._lmp = None

        self.assertTrue( hasattr( steam.libmesh, '_lmp' ) )

        with self.assertRaises( RuntimeError ):
            self.db.to_common_mesh( m_target )

    def test_trans_to_self( self ):
        """ Confirm that to_common_mesh works on a mesh already in the db
        """

        m_target   = steam.mesh.mk_cube( 1, 11 )
        self.assertFalse( hasattr( steam.libmesh, '_lmp' ) )

        target_id = self.db.to_common_mesh( m_target )

        ### Check that meshid column is uniform and matches the target_id
        self.assertTrue( (self.db.data.meshid == target_id).all() )

        ### Check each solution.mesh against the target mesh
        self.assertEqual( m_target, self.db.meshes[ target_id ] )
        for solnid in self.db.data.solnid:
            self.assertEqual( self.db.solns[solnid].mesh, m_target )

    def test_serial( self ):
        """ Run the grid conversion on one processor.
        """

        self.assertFalse( hasattr( steam.libmesh, '_lmp' ) )

        ### Define target mesh and perform conversion
        m_target   = steam.mesh.mk_cube( 1, 10 )
        target_id = self.db.to_common_mesh( m_target, nproc = 1 )

        ### Check that meshid column is uniform and matches the target_id
        self.assertTrue( (self.db.data.meshid == target_id).all() )

        ### Check each solution.mesh against the target mesh
        self.assertEqual( m_target, self.db.meshes[ target_id ] )
        for solnid in self.db.data.solnid:
            self.assertEqual( self.db.solns[solnid].mesh, m_target )

    def test_parallel( self ):
        """ Run the grid conversion on all available processors.
        """

        self.assertFalse( hasattr( steam.libmesh, '_lmp' ) )

        ### Define target mesh and perform conversion
        m_target   = steam.mesh.mk_cube( 1, 10 )

        target_id = self.db.to_common_mesh( m_target )
            
        ### Check that meshid column is uniform and matches the target_id
        self.assertTrue( (self.db.data.meshid == target_id).all() )

        ### Check each solution.mesh against the target mesh
        self.assertEqual( m_target, self.db.meshes[ target_id ] )
        for solnid in self.db.data.solnid:
            self.assertEqual( self.db.solns[solnid].mesh, m_target )

#! Run the tests
if __name__ == '__main__':
    ut.main()
