#! /bin/env python
# test_static_soln.py

import unittest as ut
import steam

#from code import interact
#interact( local = dict( globals(), **locals() ) )

class TestSubset( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Make a simple cube test case to work with.
        """

        ### Define a cube
        cls.cube = steam.mesh.mk_cube( 1, 5)
        cls.cube.get_xyz_el()
        cls.cube.static_soln = steam.solution.uniform_soln( cls.cube, 
                        node_elem_flag = 'ELEM', values = [1.0,'rand','rand'])

    def test_subset( self ):
        """Confirm that mesh.Mesh.return_subset works properly"""

        ### For whatever reason, the face components are not numbered in the 
        ### order that they're created/stored
        comp_order = [1,3,2,4,5,6]
        for i, comp in enumerate(comp_order):
            el_list = self.cube.get_comp( [comp] )
            face = self.cube.return_subset( el_list )

            face.get_xyz_el()

            ### Check the mesh that's returned
            #   xyz_el
            self.assertTrue( 
                self.cube.xyz_el.loc[face.xyz_el.index].equals( face.xyz_el ) )
            ### In Mesh.return_subset the nodes are reindexed to start from 0
            #   and increase contiguously.  Because the values in the 
            #   Mesh.conn dataframe are node IDs, slicing the full conn 
            #   dataframe will not be equal to the entire face.conn dataframe.
            #   Therefore we can only compare the indices.
            self.assertTrue( face.conn.index.equals( face.xyz_el.index ) )


            ### Check the static_soln returned with that mesh
            self.assertEqual( self.cube.static_soln.store_at,
                              face.static_soln.store_at )
            if self.cube.static_soln.store_at == 'ELEMENTS':
                self.assertTrue( 
                    self.cube.static_soln.data.loc[face.xyz_el.index].equals(
                    face.static_soln.data ) )
            elif self.cube.static_soln.store_at == 'NODES':
                raise NotImplementedError( 'Yet to develop test on node-based '
                                         + 'static_soln.' )
            else:
                raise ValueError( 'self.cube.static_soln.store_at must be '
                      + 'either "NODES" or "ELEMENTS".\nCurrent value:    '
                      + f'{self.cube.static_soln.store_at}' )

        ### Confirm that every pair of consecutively numbered elements forms 
        ### a subset that contains only 4 nodes.
        for elem in range( 0, self.cube.conn.shape[0], 2 ):
            el_list = [elem, elem+1]
            small = self.cube.return_subset( el_list )
            self.assertEqual( small.xyz_pt.shape[0], 4 )

    def test_return_comp( self ):
        """Confirm that mesh.Mesh.return_comp works properly"""

        ### For whatever reason, the face components are not numbered in the 
        ### order that they're created/stored
        comp_order = [1,3,2,4,5,6]
        for i, comp in enumerate(comp_order):
            face = self.cube.return_comp( comp_list = [comp] )

            ### Check the mesh that's returned
            face.get_xyz_el()
            self.assertTrue( 
                self.cube.xyz_el.loc[face.xyz_el.index].equals( face.xyz_el ) )

            ### Check the static_soln returned with that mesh
            self.assertEqual( self.cube.static_soln.store_at,
                              face.static_soln.store_at )
            if self.cube.static_soln.store_at == 'ELEMENTS':
                self.assertTrue( 
                    self.cube.static_soln.data.loc[face.xyz_el.index].equals(
                    face.static_soln.data ) )
            elif self.cube.static_soln.store_at == 'NODES':
                raise NotImplementedError( 'Yet to develop test on node-based '
                                         + 'static_soln.' )
            else:
                raise ValueError( 'self.cube.static_soln.store_at must be '
                      + 'either "NODES" or "ELEMENTS".\nCurrent value:    '
                      + f'{self.cube.static_soln.store_at}' )


#    @classmethod
#    def tearDownClass( cls ):
#        """Clean out temporary files."""
#
#        import os
#
#        files_to_remove = ( 'test.h5', 'node_el.pkl', 'node_no_el.pkl',
#                            'element_no_el.pkl', 'element_el.pkl' )
#        for f in files_to_remove:
#            if os.path.isfile( f ):
#                os.remove( f )

#! Run the tests
if __name__ == '__main__':
    ut.main()

