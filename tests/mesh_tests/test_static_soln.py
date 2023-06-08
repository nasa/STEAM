#! /bin/env python
# test_static_soln.py

import unittest as ut
import steam
from code import interact

class TestProject( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Make a simple cube test case to work with.
        """

        ### Cubes with solution at NODES
        #   Define a mesh that DOES NOT have xyz_el populated
        cls.node_no_el = steam.mesh.mk_cube( 1, 5)
        cls.node_no_el.static_soln = steam.solution.uniform_soln( 
                        cls.node_no_el, node_elem_flag = 'NODE', 
                        values = [1.0,'rand','rand'])
        #   Define a mesh that DOES have xyz_el populated
        cls.node_el = steam.mesh.mk_cube( 1, 5)
        cls.node_el.static_soln = steam.solution.uniform_soln( 
                        cls.node_el, node_elem_flag = 'NODE', 
                        values = [2.0,'rand','rand'])
        cls.node_el.get_xyz_el()

        #   Define a mesh that DOES NOT have xyz_el populated
        cls.element_no_el = steam.mesh.mk_cube( 1, 5)
        cls.element_no_el.static_soln = steam.solution.uniform_soln( 
                        cls.element_no_el, node_elem_flag = 'ELEMENT', 
                        values = [1.0,'rand','rand'])
        #   Define a mesh that DOES have xyz_el populated
        cls.element_el = steam.mesh.mk_cube( 1, 5)
        cls.element_el.static_soln = steam.solution.uniform_soln( 
                        cls.element_el, node_elem_flag = 'ELEMENT', 
                        values = [2.0,'rand','rand'])
        cls.element_el.get_xyz_el()

        ### Define list of files that will need to be removed
        cls.files_to_remove = [ 'node_el.pkl', 'node_no_el.pkl',
                'element_no_el.pkl', 'element_el.pkl' ]

    def test_hdf5( self ):
        """Write to an h5 file, read it, and compare.  """

        h5_file = 'test.h5'
        self.files_to_remove.append( h5_file )

        ### Write out to HDF5
        with steam.container.Container(h5_file, 'w' ) as out_cont:
            out_cont.add("/node_no_el",self.node_no_el)
            out_cont.add("/node_el",self.node_el)
            out_cont.add("/element_no_el",self.element_no_el)
            out_cont.add("/element_el",self.element_el)
        #out_cont.write()

        ### Read the HDF5
        with steam.container.Container( h5_file, 'r' ) as in_cont:
            in_cont.read()
            in_node_no_el = in_cont.obj['node_no_el']
            in_node_el = in_cont.obj['node_el']
            in_element_no_el = in_cont.obj['element_no_el']
            in_element_el = in_cont.obj['element_el']

        ### Compare node-based meshes
        #   No xyz_el
        self.assertEqual( self.node_no_el.static_soln,
                          in_node_no_el.static_soln )
        self.assertEqual( self.node_no_el, in_node_no_el )
        #   With xyz_el
        self.assertEqual( self.node_el.static_soln, in_node_el.static_soln )
        self.assertEqual( self.node_el, in_node_el )

        ### Compare element-based meshes
        #   No xyz_el
        self.assertEqual( self.element_no_el.static_soln,
                          in_element_no_el.static_soln )
        self.assertEqual( self.element_no_el, in_element_no_el )
        #   With xyz_el
        self.assertEqual( self.element_el.static_soln, 
                          in_element_el.static_soln )
        self.assertEqual( self.element_el, in_element_el )

    def test_pickle( self ):
        """Write to a pickle file, read it, and compare.  """
        
        ### Write the four meshes to pickle files
        steam.io.write_pickle( mesh = self.node_el, name = 'node_el.pkl' )
        steam.io.write_pickle( mesh = self.node_no_el, name = 'node_no_el.pkl')
        steam.io.write_pickle( mesh = self.element_el, name = 'element_el.pkl')
        steam.io.write_pickle( mesh = self.element_no_el, 
                               name = 'element_no_el.pkl' )

        ### Read and test the four pickle files
        #   Node-based mesh WITH xyz_el
        in_mesh, in_soln = steam.io.read_pickle( 'node_el.pkl' )
        self.assertTrue( in_soln is None )
        self.assertEqual( self.node_el.static_soln,
                          in_mesh.static_soln )
        self.assertEqual( self.node_el, in_mesh )

        #   Node-based mesh WITHOUT xyz_el
        in_mesh, in_soln = steam.io.read_pickle( 'node_no_el.pkl' )
        self.assertTrue( in_soln is None )
        self.assertEqual( self.node_no_el.static_soln,
                          in_mesh.static_soln )
        self.assertEqual( self.node_no_el, in_mesh )

        #   Element-based mesh WITH xyz_el
        in_mesh, in_soln = steam.io.read_pickle( 'element_el.pkl' )
        self.assertTrue( in_soln is None )
        self.assertEqual( self.element_el.static_soln,
                          in_mesh.static_soln )
        self.assertEqual( self.element_el, in_mesh )

        #   Element-based mesh WITHOUT xyz_el
        in_mesh, in_soln = steam.io.read_pickle( 'element_no_el.pkl' )
        self.assertTrue( in_soln is None )
        self.assertEqual( self.element_no_el.static_soln,
                          in_mesh.static_soln )
        self.assertEqual( self.element_no_el, in_mesh )

    def test_copy( self ):
        """Confirm that copying works properly """

        ### Copy the node-based mesh with xyz_el
        copied_node_el = steam.mesh.Mesh()
        copied_node_el.copy( self.node_el )
        #   Compare
        self.assertEqual( self.node_el.static_soln, copied_node_el.static_soln)
        self.assertEqual( self.node_el, copied_node_el )

        ### Copy the node-based mesh with no xyz_el
        copied_node_no_el = steam.mesh.Mesh()
        copied_node_no_el.copy( self.node_no_el )
        #   Compare
        self.assertEqual( self.node_no_el.static_soln, 
                          copied_node_no_el.static_soln)
        self.assertEqual( self.node_no_el, copied_node_no_el )

        ### Copy the element-based mesh with xyz_el
        copied_element_el = steam.mesh.Mesh()
        copied_element_el.copy( self.element_el )
        #   Compare
        self.assertEqual( self.element_el.static_soln, 
                          copied_element_el.static_soln)
        self.assertEqual( self.element_el, copied_element_el )

        ### Copy the element-based mesh with no xyz_el
        copied_element_no_el = steam.mesh.Mesh()
        copied_element_no_el.copy( self.element_no_el )
        #   Compare
        self.assertEqual( self.element_no_el.static_soln, 
                          copied_element_no_el.static_soln)
        self.assertEqual( self.element_no_el, copied_element_no_el )

    def test_delete_static( self ):
        """Delete a static_soln from an HDF5 file. """

        h5_file = 'test_delete.h5'
        self.files_to_remove.append( h5_file )

        ### Write out to HDF5
        with steam.container.Container( h5_file, 'w' ) as out_cont:
            out_cont.add("/node_no_el",self.node_no_el)
            out_cont.add("/node_el",self.node_el)
            out_cont.add("/element_no_el",self.element_no_el)
            out_cont.add("/element_el",self.element_el)

        ### Create another mesh that never had a static_soln 
        never_had_static = steam.mesh.mk_cube( 1, 5)
        self.assertIsNone( never_had_static.static_soln )

        ### Read the HDF5 in "r+" mode and delete all static solutions
        with steam.container.Container( h5_file, 'r+' ) as in_cont:
            in_cont.read()
            in_node_no_el = in_cont.obj['node_no_el']
            self.assertIsNotNone( in_node_no_el.static_soln )
            in_node_no_el.static_soln = None
            self.assertIsNone( in_node_no_el.static_soln )

            in_node_el = in_cont.obj['node_el']
            self.assertIsNotNone( in_node_el.static_soln )
            in_node_el.static_soln = None
            self.assertIsNone( in_node_el.static_soln )

            in_element_no_el = in_cont.obj['element_no_el']
            self.assertIsNotNone( in_element_no_el.static_soln )
            in_element_no_el.static_soln = None
            self.assertIsNone( in_element_no_el.static_soln )

            in_element_el = in_cont.obj['element_el']
            self.assertIsNotNone( in_element_el.static_soln )
            in_element_el.static_soln = None
            self.assertIsNone( in_element_el.static_soln )
    
            ### Add the mesh that never had a static_soln
            in_cont.add( "never_had", never_had_static )

        ### Read HDF5 in "r" mode and confirm static solutions aren't present
        with steam.container.Container( h5_file, 'r' ) as cont3:
            cont3.read()

            in_never = cont3.obj['never_had']
            self.assertIsNone( in_never.static_soln )

            in_node_no_el = cont3.obj['node_no_el']
            self.assertIsNone( in_node_no_el.static_soln )

            in_node_el = in_cont.obj['node_el']
            self.assertIsNone( in_node_el.static_soln )

            in_element_no_el = in_cont.obj['element_no_el']
            self.assertIsNone( in_element_no_el.static_soln )

            in_element_el = in_cont.obj['element_el']
            self.assertIsNone( in_element_el.static_soln )

    @classmethod
    def tearDownClass( cls ):
        """Clean out temporary files."""

        import os

        for f in cls.files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )

#! Run the tests
if __name__ == '__main__':
    ut.main()

