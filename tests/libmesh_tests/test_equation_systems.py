#!/usr/bin/env python
### test_equation_systems.py
###
### Test the basic functionality of the Table class.
### In order to run these tests with verbose output:
### >>> python test_equation_systems.py --verbose
### >>> python test_equation_systems.py -vb "I [Cam] believes what this does is suppresses output of tests that pass, but if it fails, it prints the output to the screen, which is a pretty nice feature!"

import unittest
import os
import steam

@unittest.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestEquationSystems( unittest.TestCase ):

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

        ### Create a LibMeshInit object
        # NOTE: this hasn't yet been required in STEAM, but it will be soon

        ### Prepare to write cube to XDR file
        self.out_xdr = "./lcube.xdr"
        self.remove_files.append( self.out_xdr )

        #! Convert cube to libmesh
        lcube = steam.libmesh.Mesh.from_python(self.cube)

        # test that we have a valid pointer
        self.assertIsNotNone(lcube())
        lcube.write( self.out_xdr )

    def test_equation_systems( self ):

        # test that we have a valid pointer
        self.assertIsNotNone(steam.libmesh._lmp)

        # create and read LibMesh Mesh object through SWIG inteface
        mesh = steam.libmesh.Mesh( self.out_xdr )

        # test that we have a valid pointer
        self.assertIsNotNone(mesh())

        #create a LibMesh EquationSystems object through SWIG inteface
        equation_systems = steam.libmesh.EquationSystems(mesh)

        # test that we have a valid pointer
        self.assertIsNotNone(equation_systems())

        # add a transient explicit system to the EquationSystems object
        equation_systems.add_transient_explicit_system("MyTransientExplicitSystem")

        # initialize EquationSystems object for use
        equation_systems.init()

        # print information
        equation_systems.print_info()

        # clean up the equation_systems
        del( equation_systems )

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
