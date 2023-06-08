#! /bin/env python

import unittest as ut
import steam
import numpy as np
import pandas as pd

def get_simple_cube(dim,l=1.0):
    """Make a cube with each side a different value."""
    mesh = steam.mesh.mk_cube( l, dim)
    mesh.get_xyz_el()
    soln = steam.solution.uniform_soln( mesh, values = [1.0])
    steam.models.scalar_to_component(soln,1,'q1',1.0)
    steam.models.scalar_to_component(soln,2,'q1',2.0)
    steam.models.scalar_to_component(soln,3,'q1',3.0)
    steam.models.scalar_to_component(soln,4,'q1',4.0)
    steam.models.scalar_to_component(soln,5,'q1',5.0)
    steam.models.scalar_to_component(soln,6,'q1',6.0)
    return (mesh,soln)

class TestProject( ut.TestCase ):

    def test_analytical( self ):
        """ This test will use a transform and do the math by hand to show that they
            are the same.  Only using 2D grid.
        """

        ### Start with three points in the source grid, their quantity is their index.
        ### The third point is a bit farther away to make sure that it isn't used when
        ### only using two points.
        s_locs = np.array([[0   ,0,0],
                           [1.01,0,0],
                           [0.5 ,1,0]])

        s_vals = np.array([[1],
                           [2],
                           [3]])


        for y_val in [0,0.25,0.5,1.0]:

        ### The targets are just a simple list of points between A and B
            t_locs = np.array([[0.0,y_val,0],
                               [0.2,y_val,0],
                               [0.4,y_val,0],
                               [0.5,y_val,0],
                               [0.6,y_val,0],
                               [0.8,y_val,0],
                               [1.01,y_val,0]])

        ### Make these numpy arrays into STEAM meshes/solutions
            s_mesh = steam.mesh.Mesh()
            s_mesh.xyz_pt = pd.DataFrame(s_locs,columns=["X","Y","Z"])
            s_mesh.to_point_cloud()

            s_soln = steam.solution.Solution()
            s_soln.mesh = s_mesh
            s_soln.data = pd.DataFrame(s_vals,columns=["Q"])

            t_mesh = steam.mesh.Mesh()
            t_mesh.xyz_pt = pd.DataFrame(t_locs,columns=["X","Y","Z"])
            t_mesh.to_point_cloud()


            for k in [1,2,3]:
                for p in [1,2,3,4]:
                    invdist = steam.interpolate.Transform(s_mesh,t_mesh)
                    print("{} :  k = {}, p = {} ...".format(y_val,k,p))
                    ### Interpolate between the first two points
                    t_vals = np.zeros((t_locs.shape[0],1))
                    for (i,pt) in enumerate(t_locs):
                        #! Fun math to get the distances to each point:
                        d_vect = s_locs - pt
                        distances = np.sqrt(np.einsum('...j,...j',d_vect,d_vect))
                        ranking = np.argsort(distances)
                        
                        closest_dist = np.zeros(k)
                        closest_vals = np.zeros(k)

                        for ki in range(k):
                            closest_dist[ki] = distances[ranking[ki]]
                            closest_vals[ki] = s_vals[   ranking[ki]]
                        ## Ok, now we are size k

                        # The last one is different since it's on the point
                        if closest_dist[0] == 0.0:
                            t_vals[i] = closest_vals[0]

                        else:
                            weights   = closest_dist**-p
                            sum_w     = np.sum(weights)
                            t_vals[i] = np.einsum('i...,i',closest_vals,weights)/sum_w

                    dist = invdist.inverse_distance(k=k,p=p)  # Only use the first two points
                    t_soln = invdist.apply(s_soln)
                    np.testing.assert_allclose(t_vals,t_soln.data.values)

                    # QUAD breaks for k > 2 since there are times when two points are in the same
                    # quadrent.  In that case, you're not going to get the same answer as the
                    # analytical since it uses all the points
                    if (k <= 2):
                        dist = invdist.inverse_distance_quad(k=k,p=p)  # Only use the first two points
                        t_soln = invdist.apply(s_soln)
                        np.testing.assert_allclose(t_vals,t_soln.data.values)
        return


    def test_exact_cube( self):
        """Interpolate from a cube to itself.  This should be trivial, but didn't always work!"""


        (mesh,soln) = get_simple_cube(5)

        target = steam.mesh.mk_cube( 1, 5 )
        target.get_xyz_el()

        for k in [1,3,5]:
            for p in [1,3,5]:
                trans =  steam.interpolate.Transform(mesh,target)
                dist  =  trans.inverse_distance(k=5,p=2)
                t_soln = trans.apply(soln)
                check =np.all(soln.data.values == t_soln.data.values)
                self.assertTrue( check )

                trans =  steam.interpolate.Transform(mesh,target)
                dist  =  trans.inverse_distance_quad(k=5,p=2)
                t_soln = trans.apply(soln)
                check =np.all(soln.data.values == t_soln.data.values)
                self.assertTrue( check )

        return


    def test_cube_sphere( self):
        """Interpolate from a cube to a sphere."""

        (mesh,soln) = get_simple_cube(15)
        target = steam.mesh.mk_sphere(.65,15)
        target.get_xyz_el()
        
        trans =  steam.interpolate.Transform(mesh,target)
        dist  =  trans.inverse_distance(k=4,p=2)
        t_soln = trans.apply(soln)
        
        # The portions of the sphere that poke out should be the values on the faces
        np.testing.assert_allclose(t_soln.data['q1'].values[target.xyz_el['X'].values < -.5],1.)
        np.testing.assert_allclose(t_soln.data['q1'].values[target.xyz_el['X'].values >  .5],2.)
        np.testing.assert_allclose(t_soln.data['q1'].values[target.xyz_el['Y'].values < -.5],3.)
        np.testing.assert_allclose(t_soln.data['q1'].values[target.xyz_el['Y'].values >  .5],4.)
        np.testing.assert_allclose(t_soln.data['q1'].values[target.xyz_el['Z'].values < -.5],5.)
        np.testing.assert_allclose(t_soln.data['q1'].values[target.xyz_el['Z'].values >  .5],6.)

        return

    def test_append_additional_k_values( self):
        """ Test that we can have components with different interpolation strategies on the same grid. """

        (mesh,soln)  = get_simple_cube(11)
        soln = steam.solution.uniform_soln( mesh, values = ['rand'])
        (target,foo) = get_simple_cube(21,l=0.6)

        #! Make the target grid crazy (diverse stencils)
        target.transform(["m0x","m0y","m0z"])
        target.transform(["r10x","r-30y","r15z"])
        target.transform(["t0.05x"])

        trans  =  steam.interpolate.Transform(mesh,target)
        test1  =  steam.interpolate.Transform(mesh,target)
        test2  =  steam.interpolate.Transform(mesh,target)
        test3  =  steam.interpolate.Transform(mesh,target)
        test4  =  steam.interpolate.Transform(mesh,target)
        test5  =  steam.interpolate.Transform(mesh,target)
        test6  =  steam.interpolate.Transform(mesh,target)

        #! Make uniform transformations, and singular ones
        dist   =  test1.inverse_distance(k=1,p=2)
        dist   =  trans.inverse_distance(k=1,p=2,target_comp=["1"])
        dist   =  test2.inverse_distance(k=2,p=2)
        dist   =  trans.inverse_distance(k=2,p=2,target_comp=["2"])
        dist   =  test3.inverse_distance(k=3,p=2)
        dist   =  trans.inverse_distance(k=3,p=2,target_comp=["3"])
        dist   =  test4.inverse_distance_quad(k=1,p=2)
        dist   =  trans.inverse_distance_quad(k=1,p=2,target_comp=["4"])
        dist   =  test5.inverse_distance_quad(k=2,p=2)
        dist   =  trans.inverse_distance_quad(k=2,p=2,target_comp=["5"])
        dist   =  test6.inverse_distance_quad(k=3,p=2)
        dist   =  trans.inverse_distance_quad(k=3,p=2,target_comp=["6"])

        soln1  =  test1.apply(soln)
        soln2  =  test2.apply(soln)
        soln3  =  test3.apply(soln)
        soln4  =  test4.apply(soln)
        soln5  =  test5.apply(soln)
        soln6  =  test6.apply(soln)
        t_soln =  trans.apply(soln)

        #! Test it
        np.testing.assert_allclose(t_soln.return_comp(["1"]).data.values,
                                    soln1.return_comp(["1"]).data.values)
        np.testing.assert_allclose(t_soln.return_comp(["2"]).data.values,
                                    soln2.return_comp(["2"]).data.values)
        np.testing.assert_allclose(t_soln.return_comp(["3"]).data.values,
                                    soln3.return_comp(["3"]).data.values)
        np.testing.assert_allclose(t_soln.return_comp(["4"]).data.values,
                                    soln4.return_comp(["4"]).data.values)
        np.testing.assert_allclose(t_soln.return_comp(["5"]).data.values,
                                    soln5.return_comp(["5"]).data.values)
        np.testing.assert_allclose(t_soln.return_comp(["6"]).data.values,
                                    soln6.return_comp(["6"]).data.values)

        #steam.io.write_vtu(name="test.vtu",soln=t_soln)
        return




#! Run the tests
if __name__ == '__main__':
    ut.main()

