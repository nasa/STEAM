#! /bin/env python

import unittest as ut
import steam
import numpy as np
import pandas as pd

def get_simple_cube(dim,l=1.0):
    """Make a cube with Pressure on the min sides."""
    mesh = steam.mesh.mk_cube( l, dim)
    mesh.get_xyz_el()
    soln = steam.solution.uniform_soln( mesh, values = [1.0])
    return (mesh,soln)

class TestAreas( ut.TestCase ):

    def test_projected_area( self ):
        """ Test projected areas for entire cube and components
        """

        # Just integrate the faces
        (mesh,soln) = get_simple_cube(11)

        # Total area should be zero
        areas = mesh.get_areas()
        np.testing.assert_allclose(areas,[0.0,0.0,0.0],atol=1e-14)

        # Check the various components (1/2 X-max/min, 3/4 -Y, etc..)
        areas = mesh.get_areas(comps=[1])
        np.testing.assert_allclose(areas, [-1.0,0.0,0.0],atol=1e-14)

        areas = mesh.get_areas(comps=[2])
        np.testing.assert_allclose(areas, [ 1.0,0.0,0.0],atol=1e-14)

        areas = mesh.get_areas(comps=[3])
        np.testing.assert_allclose(areas, [0.0,-1.0,0.0],atol=1e-14)

        areas = mesh.get_areas(comps=[4])
        np.testing.assert_allclose(areas, [0.0, 1.0,0.0],atol=1e-14)

        areas = mesh.get_areas(comps=[5])
        np.testing.assert_allclose(areas, [0.0,0.0,-1.0],atol=1e-14)

        areas = mesh.get_areas(comps=[6])
        np.testing.assert_allclose(areas, [0.0,0.0, 1.0],atol=1e-14)

        # Now try some combinations
        areas = mesh.get_areas(comps=[1,2])
        np.testing.assert_allclose(areas, [ 0.0,0.0,0.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[3,4])
        np.testing.assert_allclose(areas, [ 0.0,0.0,0.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[5,6])
        np.testing.assert_allclose(areas, [ 0.0,0.0,0.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[2,4,6])
        np.testing.assert_allclose(areas, [ 1.0, 1.0, 1.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[1,3,5])
        np.testing.assert_allclose(areas, [-1.0,-1.0,-1.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[1,2,3,5])
        np.testing.assert_allclose(areas, [ 0.0,-1.0,-1.0],atol=1e-14)

        return

    def test_wetted_area( self ):
        """ Test projected areas for entire cube and components
        """

        # Just integrate the faces
        (mesh,soln) = get_simple_cube(11)

        # Total area should be zero
        areas = mesh.get_areas(wet=True)
        np.testing.assert_allclose(areas,[2.0,2.0,2.0,6.0],atol=1e-14)

        # Check the various components (1/2 X-max/min, 3/4 -Y, etc..)
        areas = mesh.get_areas(comps=[1],wet=True)
        np.testing.assert_allclose(areas, [1.0,0.0,0.0,1.0],atol=1e-14)

        areas = mesh.get_areas(comps=[2],wet=True)
        np.testing.assert_allclose(areas, [ 1.0,0.0,0.0,1.0],atol=1e-14)

        areas = mesh.get_areas(comps=[3],wet=True)
        np.testing.assert_allclose(areas, [0.0,1.0,0.0,1.0],atol=1e-14)

        areas = mesh.get_areas(comps=[4],wet=True)
        np.testing.assert_allclose(areas, [0.0,1.0,0.0,1.0],atol=1e-14)

        areas = mesh.get_areas(comps=[5],wet=True)
        np.testing.assert_allclose(areas, [0.0,0.0,1.0,1.0],atol=1e-14)

        areas = mesh.get_areas(comps=[6],wet=True)
        np.testing.assert_allclose(areas, [0.0,0.0,1.0,1.0],atol=1e-14)

        # Now try some combinations
        areas = mesh.get_areas(comps=[1,2],wet=True)
        np.testing.assert_allclose(areas, [ 2.0,0.0,0.0,2.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[3,4],wet=True)
        np.testing.assert_allclose(areas, [ 0.0,2.0,0.0,2.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[5,6],wet=True)
        np.testing.assert_allclose(areas, [ 0.0,0.0,2.0,2.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[2,4,6],wet=True)
        np.testing.assert_allclose(areas, [ 1.0, 1.0, 1.0,3.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[1,3,5],wet=True)
        np.testing.assert_allclose(areas, [ 1.0, 1.0, 1.0,3.0],atol=1e-14)
            
        areas = mesh.get_areas(comps=[1,2,3,5],wet=True)
        np.testing.assert_allclose(areas, [ 2.0, 1.0, 1.0,4.0],atol=1e-14)

        return

#! Run the tests
if __name__ == '__main__':
    ut.main()
