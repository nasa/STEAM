#!/usr/bin/env python
### test_modified_newtonian.py
###
### Test the basic functionality of the Table class.
### In order to run these tests with verbose output:
### >>> python test_modified_newtonian.py --verbose

import unittest
import os
import steam

@unittest.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestModifiedNewtonian( unittest.TestCase ):

    @classmethod
    def setUpClass( cls ):
        cls.remove_files = []
        steam.libmesh.Init()

    def test_modified_newtonian( self ):

        # test that we have a valid pointer
        self.assertIsNotNone(steam.libmesh._lmp)

        # create and read LibMesh Mesh object through SWIG inteface
        #mesh = steam.libmesh.Mesh("../libmesh_tests/orion_cm607f.xdr")
        mesh = steam.mesh.mk_sphere( 1, 31 )
        lmesh = steam.libmesh.Mesh.from_python( mesh )

        # test that we have a valid pointer
        self.assertIsNotNone(lmesh())

        # create the freestream vector -- roughly 18 degrees alpha
        fs = [ -0.95106, 0.0, 0.30902 ]

        out_file = './modified_newtonian.plt'
        #self.remove_files.append( out_file )
        steam.physics_models.modified_newtonian.get_global_solution(
                    lmesh(),  # pointer to the mesh
                    6000.0,   # velocity m/se
                    1.0,      # density kg/m^3
                    300.0,    # T K
                    1.4,      # gamma (ratio of specific heats)
                    287.15,   # gas constant K/(kg*K)
                    fs,       # freestream vector
                    out_file ) # output filename
        
        self.assertTrue( os.path.isfile( out_file ) )

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
