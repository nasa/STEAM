#!/usr/bin/env python
### test_mpicomm_init.py
###
### In order to run these tests with verbose output:
### >>> python test_mpicomm_init.py --verbose
### >>> python test_mpicomm_init.py -vb "I [Cam] believes what this does is suppresses output of tests that pass, but if it fails, it prints the output to the screen, which is a pretty nice feature!"

import unittest
import steam

@unittest.skip( "Skipping until we know how to test MPI.  See issue #24." )
class TestInit( unittest.TestCase ):

    def setUp( self ):
        ### Create a mesh on the Python side with which to test LibMesh
        self.mesh = steam.mesh.mk_cube( 1, 5 )

        ### Create a new MPI_Comm the pass to LibMeshInit

    def test_mpicomm_init( self ):

        ### Test that we don't have a valid LibMeshInit object
        self.assertRaises( NameError, steam.libmesh.Mesh )
        self.assertFalse( hasattr(steam.libmesh, "_lmp") )
    
        ### Create an init
        init = steam.libmesh.Init()

        ### Test that we now have a valid LibMeshInit object
        self.assertTrue( hasattr(steam.libmesh, "_lmp") )
        lmesh = steam.libmesh.Mesh.from_python( self.mesh )
        self.assertIsNotNone(steam.libmesh._lmp)
        self.assertIs( init, steam.libmesh._lmp )
        self.assertIs( init, lmesh.parent )
    

### Run the tests
if __name__ == '__main__':
    unittest.main()
