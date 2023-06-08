#! /bin/env python

import unittest as ut
import steam
import numpy as np
import pandas as pd

def get_simple_cube():
    """Make a cube with each side a different value."""
    mesh = steam.mesh.mk_cube( 1, 21 )
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

    def test_source_subset( self):
        """Interpolate from a component on a cube to another cube"""

        (s_mesh,s_soln) = get_simple_cube()

        t_mesh = steam.mesh.mk_cube( 1, 11)
        t_mesh.get_xyz_el()

        #! OK, interpolate from only one component of the cube
        trans  =  steam.interpolate.Transform(s_mesh,t_mesh)
        dist   = trans.inverse_distance(k=3,source_comp=['1'])
        t_soln = trans.apply(s_soln)
            #! So, everything should be 1.0 since it's from side 1.0
        np.testing.assert_allclose(t_soln.data.values,1.0)

        #! Interpolate from more than one component
            #! Side 1 is xmin, Side 2 is xmax, so they should each
            #! collor one half of the cube
        trans  =  steam.interpolate.Transform(s_mesh,t_mesh)
        dist   = trans.inverse_distance(k=3,source_comp=['1','2'])
        t_soln = trans.apply(s_soln)
            #! Negative side is from component 2
        np.testing.assert_allclose(t_soln.data.values[t_mesh.xyz_el['X'].values < 0,0],1.0)
            #! Positive side is from component 2
        np.testing.assert_allclose(t_soln.data.values[t_mesh.xyz_el['X'].values > 0,0],2.0)

        #! Use quad interpolate and check to make sure that we're...?
        return

    def test_target_subset( self):
        """Interpolate from a cube to a component on a cube"""

        (s_mesh,s_soln) = get_simple_cube()
        
        t_mesh = steam.mesh.mk_cube( 1, 11)
        t_mesh.get_xyz_el()
        
        #! OK, interpolate from only one component of the cube
        trans  =  steam.interpolate.Transform(s_mesh,t_mesh)
        dist   = trans.inverse_distance(k=3)  # Default to all-to-all
        dist   = trans.inverse_distance(k=3,source_comp=['6'],target_comp=['1']) # Component 1 is different
        t_soln = trans.apply(s_soln)
            #! So, everything should be the same, but comp 1 is replaced
        np.testing.assert_allclose(t_soln.return_comp(['1']).data.values,6.0)   # The changed one
        np.testing.assert_allclose(t_soln.return_comp(['2']).data.values,2.0)
        np.testing.assert_allclose(t_soln.return_comp(['3']).data.values,3.0)
        np.testing.assert_allclose(t_soln.return_comp(['4']).data.values,4.0)
        np.testing.assert_allclose(t_soln.return_comp(['5']).data.values,5.0)
        np.testing.assert_allclose(t_soln.return_comp(['6']).data.values,6.0)

        #! Now do a lot of changes, swap things around
        dist   = trans.inverse_distance(k=3,source_comp=['1'],target_comp=['2'])
        dist   = trans.inverse_distance(k=3,source_comp=['2'],target_comp=['1'])
        dist   = trans.inverse_distance(k=3,source_comp=['3'],target_comp=['4'])
        dist   = trans.inverse_distance(k=3,source_comp=['4'],target_comp=['3'])
        dist   = trans.inverse_distance(k=3,source_comp=['5'],target_comp=['6'])
        dist   = trans.inverse_distance(k=3,source_comp=['6'],target_comp=['5'])
        t_soln = trans.apply(s_soln)

        np.testing.assert_allclose(t_soln.return_comp(['1']).data.values,2.0)
        np.testing.assert_allclose(t_soln.return_comp(['2']).data.values,1.0)
        np.testing.assert_allclose(t_soln.return_comp(['3']).data.values,4.0)
        np.testing.assert_allclose(t_soln.return_comp(['4']).data.values,3.0)
        np.testing.assert_allclose(t_soln.return_comp(['5']).data.values,6.0)
        np.testing.assert_allclose(t_soln.return_comp(['6']).data.values,5.0)

        return

    def test_empty_target( self):
        """Make sure we get an error if we try to interpolate using an incomplete transformation.
        """

        (s_mesh,s_soln) = get_simple_cube()
        
        t_mesh = steam.mesh.mk_cube( 1, 11)
        t_mesh.get_xyz_el()
        
        #! Interpolate for only one component on target
        trans  =  steam.interpolate.Transform(s_mesh,t_mesh)
        dist   = trans.inverse_distance(target_comp=['1']) 
        
        self.assertRaises(Exception,trans.apply,s_soln)

        return

#! Run the tests
if __name__ == '__main__':
    ut.main()


#! Run the tests
if __name__ == '__main__':
    ut.main()

