#!/usr/bin/env python
### test_mesh.py
###
### Test the basic functionality of the Table class.
### In order to run these tests with verbose output:
### >>> python test_mesh.py --verbose

import unittest
import steam
import os
import numpy as np

@unittest.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestMesh( unittest.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """ Class called only once, before all tests.
        """
        cls.remove_files = []

        ### Create a LibMeshInit object
        cls.init = steam.libmesh.Init()

    def setUp( self ):
        """ Prepare everything we need for these tests.  We'll use the same
        cube mesh for each test.
        """

        #! Make the cube
        self.cube = steam.mesh.mk_cube(5.0,31)

        ### Save a few parameters about it
        self.cube.get_xyz_el()
        self.n_pts =  self.cube.xyz_pt.shape[0]
        self.n_elem = self.cube.xyz_el.shape[0]
        self.out_xdr = "./lcube.xdr"
        self.remove_files.append( self.out_xdr )

    #@unittest.skip( 'Skip during rewrite' )
    def test_read_write( self ):
        """ Test XDR reading and writing.
        This test depends on the same functionality and therefore will 
        fail if test_python_to_libmesh fails.
        """

        # test that we have a valid pointer
        self.assertIsNotNone(steam.libmesh._lmp)

        #! Convert cube to libmesh
        lcube = steam.libmesh.Mesh.from_python(self.cube)

        # test that we have a valid pointer
        self.assertIsNotNone(lcube())

        # write the mesh in EXODUS format
        lcube.write( self.out_xdr )

        # read and read the grid we just wrote in EXODUS II format
        lcube2 = steam.libmesh.Mesh( self.out_xdr )

        # test that we have a valid pointer
        self.assertIsNotNone(lcube2())

        # Convert lcube2 back to python
        cube2 = steam.mesh.Mesh.from_libmesh(lcube2)
        cube2.get_xyz_el()

        ### Compare the old and new cube meshes
        #   NOTE: reading the XDR file creates a mesh that is different from
        #         the original.  We must therefore test to confirm that they
        #         are functionally equivalent, which can't be done
        #         using steam.mesh.__eq__
        # Confirm that all points are in the same places
        np.testing.assert_equal(
                self.cube.xyz_pt.sort_values(by=['X','Y','Z']).values,
                    cube2.xyz_pt.sort_values(by=['X','Y','Z']).values )
        # Confirm that all element centroids are in the same places
        np.testing.assert_equal(
                self.cube.xyz_el.sort_values(by=['X','Y','Z']).values,
                    cube2.xyz_el.sort_values(by=['X','Y','Z']).values )

        # do some superficial tests on the meshes to make sure we they are the same and they are what we expect
        self.assertEqual(self.n_pts,   5402)
        self.assertEqual(self.n_elem, 10800)
        self.assertEqual(cube2.xyz_pt.shape[0],  5402)
        self.assertEqual(cube2.xyz_el.shape[0], 10800)

        # delete the EXODUS II file
        os.remove( self.out_xdr )

    def test_python_to_libmesh(self):
        """ Test the python-to-libmesh converion.  """

        #! Convert it to libmesh
        lcube = steam.libmesh.Mesh.from_python(self.cube)

        #! Convert it back
        cube2 = steam.mesh.Mesh.from_libmesh(lcube)

        #! Compare everything
        self.assertEqual( self.cube, cube2 )

        ### NOTE: Everything below is compared in the line above this one,
        #         but there's not much reason to get rid of it.

        # Points are the same:
        # - first  all() collapses rows
        # - second all() collapses columns
        all_points = ((self.cube.xyz_pt == cube2.xyz_pt).all()).all()
        self.assertTrue(all_points,"Points do not match")
        # Connectivity is  the same:
        all_conn   = ((self.cube.conn   == cube2.conn  ).all()).all()
        self.assertTrue(all_conn  ,"Connectivity does not match")
        # Components   are the same:
        all_comp   = ((self.cube.comp   == cube2.comp  ).all()).all()
        self.assertTrue(all_comp  ,"Components do not match")

    def test_mesh_normals( self ):

        # test that we have a valid pointer
        self.assertIsNotNone(steam.libmesh._lmp)

        # create and read LibMesh Mesh object through SWIG inteface
        mesh = steam.libmesh.Mesh.from_python( self.cube )

        # test that we have a valid pointer
        self.assertIsNotNone(mesh())

        # calculate normals
        mesh.populate_node_normal();

        # write the normals
        node_normal_file = "./node_normals.dat"
        mesh.write_node_normals( node_normal_file )
        self.remove_files.append( node_normal_file )
        self.assertTrue( os.path.isfile( node_normal_file ) )

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """

        for f in cls.remove_files:
            if os.path.isfile( f ):
                os.remove( f )

### Run the tests
if __name__ == '__main__':
    unittest.main()
