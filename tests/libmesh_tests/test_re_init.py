#!/usr/bin/env python
### test_re_init.py
###
### In order to run these tests with verbose output:
### >>> python test_re_init.py --verbose
### >>> python test_re_init.py -vb "I [Cam] believes what this does is suppresses output of tests that pass, but if it fails, it prints the output to the screen, which is a pretty nice feature!"

import unittest
import steam

@unittest.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestInit( unittest.TestCase ):

    def test_re_init( self ):
        """ Test that we can't reinitialize LibMesh once it's already been done.
        """

        ### Test that we don't have a valid LibMeshInit object
        if hasattr(steam.libmesh, "_lmp"):
            ### Some other unit test has already initialized LibMesh
            #   Attempt to create another Init -- should raise warning
            self.assertWarns( RuntimeWarning, steam.libmesh.Init )


            init = steam.libmesh.Init()
            self.assertEqual( init.p, steam.libmesh._lmp.p )
        else:
            ### This test is run on its own; no init exists
            self.assertRaises( RuntimeError, steam.libmesh.Mesh )
    
            ### Create an Init
            init = steam.libmesh.Init()
            self.assertIsNotNone(steam.libmesh._lmp)

            self.assertIs( init, steam.libmesh._lmp )

        ### Attempt to create another Init -- should raise warning
        self.assertWarns( RuntimeWarning, steam.libmesh.Init )
        init2 = steam.libmesh.Init()
        self.assertEqual( init2.p, steam.libmesh._lmp.p )

### Run the tests
if __name__ == '__main__':
    unittest.main()
