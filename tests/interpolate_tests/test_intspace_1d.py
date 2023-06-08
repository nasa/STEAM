#! /bin/env python

import unittest as ut
import numpy as np
import steam

class TestIntSpace1D( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Create temporary files. This method gets called once, at the very
        beginning of testing.
        """
        pass

    def setUp( self ):
        """ Set-up the points and create a space.

        This method gets called before each test_* method in order to prepare the proper variables/data.
        """ 

        # Make a simple database in X
        # Variables at X, linear, log-linear, log10-linear, and square-linear

        #! Columns:
        #! - 0: X   1: 2*X    2: log(X)    3: log10(X)    4: (3*X)**2.0
        self.points = np.ndarray((4,5))
        self.points[0:4,0] = [1,10,100,1000]
        self.points[0:4,1] = self.points[0:4,0]*2.0
        self.points[0:4,2] = np.log(self.points[0:4,0])
        self.points[0:4,3] = np.log10(self.points[0:4,0])
        self.points[0:4,4] = (self.points[0:4,0]*3.0)**2.0

        #! Random order.  This bit us once.
        np.random.shuffle(self.points)

        #! Number of test points to generate for each test
        self.ntestp        = 10



    def test_SimpleInterp( self ):
        """ Just interpolate inside of the space and check the weights/indicies.  """

        #! Columns:
        #! - 0: X   1: 2*X    2: log(X)    3: log10(X)    4: (3*X)**2.0
        for i,opts in enumerate([dict(),{'X':'loge'},{'X':'log10'},{'X':'**2'}]):

            space = steam.interpolate.IntSpace(indep=['X'])
            space.build_space(self.points[:,0],opts=opts)

            #! Test output at random points
            test_points = np.random.random(self.ntestp)*1000.0

            #! Get the weights and then get the dep variables
            (idx,wei)  = space.get_weights(test_points)
            out_points = np.einsum('ijk,ij->i...k',self.points[np.array(idx)],wei)

            import sys
            function_name = sys._getframe().f_code.co_name

            if i == 0:
                print(" - Testing linear in {}".format(function_name))
                np.testing.assert_allclose(test_points,out_points[:,0],atol=1e-14)
                np.testing.assert_allclose(test_points*2.0,out_points[:,1],atol=1e-14)
                print(" - Good !")
                
            if i == 1:
                print(" - Testing log-e  in {}".format(function_name))
                np.testing.assert_allclose(np.log(test_points),out_points[:,2],atol=1e-14)
                print(" - Good !")
                
            if i == 2:
                print(" - Testing log-10 in {}".format(function_name))
                np.testing.assert_allclose(np.log10(test_points),out_points[:,3],atol=1e-14)
                print(" - Good !")
                
            if i == 3:
                print(" - Testing square in {}".format(function_name))
                np.testing.assert_allclose((3*test_points)**2.0,out_points[:,4],atol=1e-14)
                print(" - Good !")
                

    def test_SubsetInterp( self ):
        """ Just interpolate inside subset of the space and check the weights/indicies.  """

        sub = [0,1,3]

        #! Columns:
        #! - 0: X   1: 2*X    2: log(X)    3: log10(X)    4: (3*X)**2.0
        for i,opts in enumerate([dict(),{'X':'loge'},{'X':'log10'},{'X':'**2'}]):

            space = steam.interpolate.IntSpace(indep=['X'])
            space.build_space(self.points[sub,0],out_index=sub,opts=opts,scale='auto')

            #! Test output at random points - but make sure they are bounded by the database
            db_min = np.min(self.points[sub,0])
            db_max = np.max(self.points[sub,0])
            db_range = db_max - db_min
            test_points = np.random.random(self.ntestp)*db_range+db_min

            #! Get the weights and then get the dep variables
            (idx,wei)  = space.get_weights(test_points)
            out_points = np.einsum('ijk,ij->i...k',self.points[np.array(idx)],wei)

            import sys
            function_name = sys._getframe().f_code.co_name

            if i == 0:
                print(" - Testing linear in {}".format(function_name))
                np.testing.assert_allclose(test_points,out_points[:,0],atol=1e-14)
                np.testing.assert_allclose(test_points*2.0,out_points[:,1],atol=1e-14)
                print(" - Good !")
                
            if i == 1:
                print(" - Testing log-e  in {}".format(function_name))
                np.testing.assert_allclose(np.log(test_points),out_points[:,2],atol=1e-14)
                print(" - Good !")
                
            if i == 2:
                print(" - Testing log-10 in {}".format(function_name))
                np.testing.assert_allclose(np.log10(test_points),out_points[:,3],atol=1e-14)
                print(" - Good !")
                
            if i == 3:
                print(" - Testing square in {}".format(function_name))
                np.testing.assert_allclose((3*test_points)**2.0,out_points[:,4],atol=1e-14)
                print(" - Good !")
                


    def test_ReadWrite( self):
        """ Check after a read/write. """

        #! Columns:
        #! - 0: X   1: 2*X    2: log(X)    3: log10(X)    4: (3*X)**2.0
        for i,opts in enumerate([dict(),{'X':'loge'},{'X':'log10'},{'X':'**2'}]):

            space = steam.interpolate.IntSpace(indep=['X'])
            space.build_space(self.points[:,0],opts=opts,scale='auto')

            with steam.container.Container("test_is.h5","w") as cont:
                cont.add("space",space)

            with steam.container.Container("test_is.h5") as cont2:
                cont2.read()

            dspace = cont2.obj['space']

            # Just check that they self-report the same
            self.assertEqual(str(space),str(dspace))

            #! Test output at random points
            test_points = np.random.random(self.ntestp)*1000.0

            #! Get the weights and then get the dep variables
            (midx,mwei)  =  space.get_weights(test_points)
            ( idx, wei)  = dspace.get_weights(test_points)
            np.testing.assert_allclose(idx,midx)
            np.testing.assert_allclose(wei,mwei)
            out_points = np.einsum('ijk,ij->i...k',self.points[np.array(idx)],wei)

            import sys
            function_name = sys._getframe().f_code.co_name

            if i == 0:
                print(" - Testing linear in {}".format(function_name))
                np.testing.assert_allclose(test_points,out_points[:,0],atol=1e-14)
                np.testing.assert_allclose(test_points*2.0,out_points[:,1],atol=1e-14)
                print(" - Good !")
                
            if i == 1:
                print(" - Testing log-e  in {}".format(function_name))
                np.testing.assert_allclose(np.log(test_points),out_points[:,2],atol=1e-14)
                print(" - Good !")
                
            if i == 2:
                print(" - Testing log-10 in {}".format(function_name))
                np.testing.assert_allclose(np.log10(test_points),out_points[:,3],atol=1e-14)
                print(" - Good !")
                
            if i == 3:
                print(" - Testing square in {}".format(function_name))
                np.testing.assert_allclose((3*test_points)**2.0,out_points[:,4],atol=1e-14)
                print(" - Good !")
                


    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """
        import os

        files_to_remove = ('test_is.h5')
        for f in files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )




#! Run the tests
if __name__ == '__main__':
    ut.main()

