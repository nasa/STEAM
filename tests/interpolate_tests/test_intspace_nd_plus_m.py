#! /bin/env python

import unittest as ut
import numpy as np
import steam
import math

class TestIntSpaceNDpM( ut.TestCase ):
    """ This is going to test the ND+M concept for interpolation.

    This assumes that the user has M 'square' parameters at which there N
    'non-square' parameters.  Interpolation finds the bounding ND databases,
    interpolates inside of them and then interpolates between them according 
    to the values in the M square parameters.
    """

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

        # Make a simple database in L,M,N, and X,Y,Z.  This will be a 3D+3 interpolation
        # where we have square coverage in L, M, and N but random coverage in X, Y, and Z.
        # Coverage in L,M,N will vary, but X,Y,Z will be a unit square
        lvals = [1,10,100,1000]
        mvals = [1,2,3,4,5]
        nvals = [-10,-5,0,5,10]

        npt_xyz = 10 # Number of randomly placed points in each XYZ sub-database

        db_points = np.array([]).reshape(0,6)
        for l in lvals:
            for m in mvals:
                for n in nvals:

                    ndim = 3 # X, Y, Z
                    bounds = np.zeros((ndim,2))
                    bounds[:,0] = 0.1
                    bounds[:,1] = 1.0

                    #! Fix the random number seed
                    np.random.seed(0)

                    from itertools import product  # This is super cool

                    #! First make the database points by using the bounds
                    for point in product(*bounds):
                        db_points = np.vstack((db_points,[l,m,n,*point]))

                    #! Now add random interior points
                    for i in range(npt_xyz):
                        point = np.random.rand(3) * 0.9 + 0.1
                        db_points = np.vstack((db_points,[l,m,n,*point]))

        # OK, Now I have my database of independant params
        # Make a set of dependent parameters with a lot of options  all are the sum of:
        #! - 0: linear   
        #!   1: loge(L),log10(X)
        #!   2: M**2,X**3,log10(Z)
        #!   3: log10(L,M) + N
        #!   4: loge(L,M) + N +(MD)**3
        dep_points = np.zeros((db_points.shape[0],5))
        dep_points[:,0] = db_points.sum(axis=1)
        dep_points[:,1] = np.log(db_points[:,0]) + np.log10(db_points[:,3])
        dep_points[:,2] = db_points[:,1]**2.0+db_points[:,3]**3.0+np.log10(db_points[:,5])
        dep_points[:,3] = np.log10(db_points[:,:2]).sum(axis=1) + db_points[:,2]
        dep_points[:,4] = np.log(db_points[:,:2]).sum(axis=1) + db_points[:,2] + (db_points[:,3:6]**3.0).sum(axis=1) 

        self.indep   = ["L","M","N","X","Y","Z"]
        self.idep_pt = db_points
        self.dep_pt  = dep_points

        #! Number of test points to generate for each test
        self.ntestp        = 10


    def test_SimpleInterp( self ):
        """ Just interpolate inside of the space and check the weights/indicies.  """

        #! Columns:
        for i,opts in enumerate([
                                 dict(),
                                 {'L':'loge','X':'log10'},
                                 {'Z':'log10','M':'**2','X':'**3'},
                                 {'M':'log10','L':'log10'},
                                 {'M':'loge','L':'loge','X':'**3','Y':'**3','Z':'**3'}
                                 ]):

            space = steam.interpolate.IntSpace(indep=self.indep,subvars=['L','M','N'])
            space.build_space(self.idep_pt,opts=opts)

            #! Test output at random points
            test_points = np.random.random((self.ntestp,6))
            test_points[:,0]   = test_points[:,0]   *999.0 + 1.0
            test_points[:,1]   = test_points[:,1]   *4.0   + 1.0
            test_points[:,2]   = test_points[:,2]   *20    - 10.0
            test_points[:,3:6] = test_points[:,3:6] * 0.9  + 0.1

            #! Get the weights and then get the dep variables
            (idx,wei)  = space.get_weights(test_points)
            out_points = np.einsum('ijk,ij->i...k',self.dep_pt[np.array(idx)],wei)
            #import code
            #code.interact(local=locals())

            import sys
            function_name = sys._getframe().f_code.co_name

            if i == 0:
                print(" - Testing case {} in {}".format(i,function_name))
                output = test_points.sum(axis=1)
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 1:
                print(" - Testing case {} in {}".format(i,function_name))
                output = np.log(test_points[:,0]) + np.log10(test_points[:,3])
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 2:
                print(" - Testing case {} in {}".format(i,function_name))
                output = test_points[:,1]**2.0+test_points[:,3]**3.0+np.log10(test_points[:,5])
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 3:
                print(" - Testing case {} in {}".format(i,function_name))
                output = np.log10(test_points[:,:2]).sum(axis=1) + test_points[:,2]
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 4:
                print(" - Testing case {} in {}".format(i,function_name))
                output = np.log(test_points[:,:2]).sum(axis=1) + test_points[:,2] + (test_points[:,3:6]**3.0).sum(axis=1) 
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
     

    def test_ReadWrite( self):
        """ Check after a read/write. 
        
        We use scale=auto here because it is important to make sure that these factors are
        invariant after a save.  If not, then it won't survive a write to disk."""

        #! Columns:
        for i,opts in enumerate([
                                 dict(),
                                 {'L':'loge','X':'log10'},
                                 {'Z':'log10','M':'**2','X':'**3'},
                                 {'M':'log10','L':'log10'},
                                 {'M':'loge','L':'loge','X':'**3','Y':'**3','Z':'**3'}
                                 ]):

            space = steam.interpolate.IntSpace(indep=self.indep,subvars=['L','M','N'])
            space.build_space(self.idep_pt,opts=opts,scale='auto')

            with steam.container.Container("test_is.h5","w") as cont:
                cont.add("space",space)

            with steam.container.Container("test_is.h5") as cont2:
                cont2.read()

            dspace = cont2.obj['space']

            # Just check that they self-report the same
            self.assertEqual(str(space),str(dspace))
 
            #! Test output at random points
            test_points = np.random.random((self.ntestp,6))
            test_points[:,0]   = test_points[:,0]   *999.0 + 1.0
            test_points[:,1]   = test_points[:,1]   *4.0   + 1.0
            test_points[:,2]   = test_points[:,2]   *20    - 10.0
            test_points[:,3:6] = test_points[:,3:6] * 0.9  + 0.1

            #! Get the weights and then get the dep variables
            (midx,mwei)  =  space.get_weights(test_points)
            (idx,wei)    = dspace.get_weights(test_points)
            np.testing.assert_allclose(midx,idx)
            np.testing.assert_allclose(mwei,wei)
            out_points = np.einsum('ijk,ij->i...k',self.dep_pt[np.array(idx)],wei)
            #import code
            #code.interact(local=locals())

            import sys
            function_name = sys._getframe().f_code.co_name

            if i == 0:
                print(" - Testing case {} in {}".format(i,function_name))
                output = test_points.sum(axis=1)
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 1:
                print(" - Testing case {} in {}".format(i,function_name))
                output = np.log(test_points[:,0]) + np.log10(test_points[:,3])
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 2:
                print(" - Testing case {} in {}".format(i,function_name))
                output = test_points[:,1]**2.0+test_points[:,3]**3.0+np.log10(test_points[:,5])
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 3:
                print(" - Testing case {} in {}".format(i,function_name))
                output = np.log10(test_points[:,:2]).sum(axis=1) + test_points[:,2]
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
                print(" - Good !")
                
            if i == 4:
                print(" - Testing case {} in {}".format(i,function_name))
                output = np.log(test_points[:,:2]).sum(axis=1) + test_points[:,2] + (test_points[:,3:6]**3.0).sum(axis=1) 
                np.testing.assert_allclose(output,out_points[:,i],atol=1e-14)
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

