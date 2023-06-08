#! /bin/env python

import unittest as ut
import steam

# Does it depend on LibMesh?
@ut.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestProject( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Create temporary files. This method gets called once, at the very
        beginning of testing.
        """
        pass

    def setUp( self ):
        """This method gets called before each test_* method in order to 
        prepare the proper variables/data.
        """ 
        pass

    def test_NAME( self ):
        """Run your tests.  You can have as many test_* methods, each testing
        various aspects of your code.
        """
        pass

    def tearDown( self ):
        """This method gets called after each test_* method in order to 
        clean up for the next test.
        """ 
        pass

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """
        pass
        #import os

        #os.remove('file')
        #import shutil
        #shutil.rmtree('dir')

        #files_to_remove = ( 'test.dat', 'test.e', 'original_soln.triq', 
        #                    'new_soln.triq', 'cube.tri' )
        #for f in files_to_remove:
        #    if os.path.isfile( f ):
        #        os.remove( f )


#! Run the tests
if __name__ == '__main__':
    ut.main()

