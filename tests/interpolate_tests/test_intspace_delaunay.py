#! /bin/env python

import unittest as ut
import numpy as np
import steam

class TestIntSpaceTri( ut.TestCase ):

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

        # Make a simple database in X,Y,Z
        # - Slanted range in X/Y.  Totally populated between X (2-100), Y(0-3)
        # - Duplicated coverage in Z (1-2).
        self.points = np.ndarray((8,3))
        self.points[0:4,0] = [1,100,1,100]
        self.points[4:8,0] = [1,100,2,200]
        self.points[0:4,1] = [0,0,3,3]
        self.points[4:8,1] = [0,0,3,3]
        self.points[0:4,2] = [1,1,1,1]
        self.points[4:8,2] = [2,2,2,2]

    def test_SimpleInterp( self ):
        """ Just interpolate inside of the space and check the weights/indicies.  """

        space = steam.interpolate.IntSpace(['X','Y','Z'])
        space.build_space(self.points)

        #! Check to make sure that I recover at the points
        for (i,p) in enumerate(self.points):
            (idx,wei) = space.get_weights(np.array([p]))
            # The weight should be 1.0 for the index in question so let's use it as a mask
            intwei = int(np.sum(idx[0]*wei[0]))
            self.assertEqual(i,intwei)

        #! Check some test points.  Should be half-way between things with linear scaling
        #! These may break if we're not on an edge.
        for w1i in range(32):
            for w2i in range(32):
                w1 = w1i/100.0
                w2 = w2i/100.0
                w3 = 1.0-w1-w2
                wei_in = np.array([w1,w2,w3])
                idx_in = np.random.randint(0,8,3)

                p = np.einsum('ij,i',self.points[idx_in],wei_in)
                (idx,wei) = space.get_weights(np.array([p]))
                p_out = np.einsum('ij,i',self.points[idx[0]],wei[0])

                np.testing.assert_allclose(p,p_out,atol=1e-14)


    def test_SubsetInterp( self ):
        """ Just interpolate inside subset of the space and check the weights/indicies.  """

        space = steam.interpolate.IntSpace(['X','Z'])
        sub   = np.array([0,1,4,5],dtype="int")
        subp  = np.ndarray((4,2))
        subp[:,0] = self.points[sub,0]
        subp[:,1] = self.points[sub,2]

        space.build_space(subp,sub)

        #! Check to make sure that I recover at the points
        for (i,p) in enumerate(subp):
            (idx,wei) = space.get_weights(np.array([p]))
            # The weight should be 1.0 for the index in question so let's use it as a mask
            intwei = int(np.sum(idx[0]*wei[0]))
            self.assertEqual(sub[i],intwei)

        #! More tests here on weights?  I think that if SimpleInterp works, then this will, too

    def test_Operations( self ):
        """ Test the built-in operations for manipulating independant parameters.  """

        space = steam.interpolate.IntSpace(['X','Y','Z'])
        opts  = {"X":"log10","Y":"**2","Z":"loge"}
        space.build_space(self.points,opts=opts)

        #! Check to make sure that I recover at the points
        for (i,p) in enumerate(self.points):
            (idx,wei) = space.get_weights(np.array([p]))
            # The weight should be 1.0 for the index in question so let's use it as a mask
            intwei = int(np.sum(idx[0]*wei[0]))
            self.assertEqual(i,intwei)

        #! Check that log10 scaling puts 10 half way between 1 and 100
        (idx,wei) = space.get_weights(np.array([[10.0,0.0,1.0]]))
        np.testing.assert_equal([0,1],sorted(idx[0][np.isclose(wei[0],0.5)]))

        #! Check that **2   works
        (idx,wei) = space.get_weights(np.array([[1.0,4.5**0.5,1.0]]))
        np.testing.assert_equal([0,2],sorted(idx[0][np.isclose(wei[0],0.5)]))

        #! Check that loge  works
        (idx,wei) = space.get_weights(np.array([[1.0,0.0,2.0**.5]]))
        np.testing.assert_equal([0,4],sorted(idx[0][np.isclose(wei[0],0.5)]))

    def test_Scaling( self ):
        """ Test the built-in operations for scaling independant parameters.  """

        for var in ['X','Y','Z']:

            space = steam.interpolate.IntSpace(['X','Y','Z'])
            opts  = {"X":"log10","Y":"**2","Z":"loge"}
            space.build_space(self.points,opts=opts,scale={var:"auto"})

            #! Check to make sure that I recover at the points
            for (i,p) in enumerate(self.points):
                (idx,wei) = space.get_weights(np.array([p]))
                # The weight should be 1.0 for the index in question so let's use it as a mask
                intwei = int(np.sum(idx[0]*wei[0]))
                self.assertEqual(i,intwei)

            #! Check that log10 scaling puts 10 half way between 1 and 100
            (idx,wei) = space.get_weights(np.array([[10.0,0.0,1.0]]))
            np.testing.assert_equal([0,1],sorted(idx[0][np.isclose(wei[0],0.5)]))

            #! Check that **2   works
            (idx,wei) = space.get_weights(np.array([[1.0,4.5**0.5,1.0]]))
            np.testing.assert_equal([0,2],sorted(idx[0][np.isclose(wei[0],0.5)]))

            #! Check that loge  works
            (idx,wei) = space.get_weights(np.array([[1.0,0.0,2.0**.5]]))
            np.testing.assert_equal([0,4],sorted(idx[0][np.isclose(wei[0],0.5)]))


    def test_ReadWrite( self):
        """ Turn on all that you can and check before/after a read/write. """


        space = steam.interpolate.IntSpace(['X','Y','Z'])
        opts  = {"X":"log10","Y":"**2","Z":"loge"}
        space.build_space(self.points,scale='auto',opts=opts)


        with steam.container.Container("test_is.h5","w") as cont:
            cont.add("space",space)

        with steam.container.Container("test_is.h5") as cont2:
            cont2.read()

        dspace = cont2.obj['space']

        # Just check that they self-report the same
        self.assertEqual(str(space),str(dspace))

        #! Check to make sure that I recover at the points
        for (i,p) in enumerate(self.points):
            ( idx, wei) =  space.get_weights(np.array([p]))
            (didx,dwei) = dspace.get_weights(np.array([p]))
            np.testing.assert_allclose(idx,didx)
            np.testing.assert_allclose(wei,dwei)

            # The weight should be 1.0 for the index in question so let's use it as a mask
            intwei = int(np.sum(idx[0]*wei[0]))
            self.assertEqual(i,intwei)

            # The weight should be 1.0 for the index in question so let's use it as a mask
            intwei = int(np.sum(didx[0]*dwei[0]))
            self.assertEqual(i,intwei)

        #! Check some test points.  Should be half-way between things with linear scaling
        #! These may break if we're not on an edge.
        for w1i in range(32):
            for w2i in range(32):
                w1 = w1i/100.0
                w2 = w2i/100.0
                w3 = 1.0-w1-w2
                wei_in = np.array([w1,w2,w3])
                idx_in = np.random.randint(0,8,3)

                p = np.einsum('ij,i',self.points[idx_in],wei_in)

                ( idx, wei) =  space.get_weights(np.array([p]))
                (didx,dwei) = dspace.get_weights(np.array([p]))
                np.testing.assert_allclose(idx,didx)
                np.testing.assert_allclose(wei,dwei)


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

