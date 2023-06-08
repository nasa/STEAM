 #! /bin/env python

import unittest as ut
import numpy as np
import steam

class TestProject( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Create temporary files. This method gets called once, at the very
        beginning of testing.
        """
        pass

    def setUp( self ):
        """ Create the necessary half-body and full-body meshes and solutions.
        """ 
        
        ### The "half" mesh will be offset by an arbitrary amount along the 
        #   x-axis and then sliced along the y-z plane.
        self.half_mesh = steam.mesh.mk_sphere( 5, 51 )
        self.offset = 2.555
        self.half_mesh.transform( [f't{offset}x'] )
        #interact( local = dict( globals(), **locals() ) )
        self.half_mesh.slice( [[(0,0,0), (0,1,0), (0,0,1)]], delete=[(-1,1)] )
        self.half_mesh.get_xyz_el()
        
        ### Create the solution for the half mesh
        self.half_soln = steam.solution.uniform_soln( self.half_mesh )
        self.half_soln.data['x10'] = self.half_mesh.xyz_el.X * 10

    @ut.skip( 'NOT YET COMPLETE' )
    def test_mirror( self ):
        """ """
        ### Make the cubes
        orig = steam.mesh.mk_cube( 1, 5)
        orig.get_xyz_el()
        cube = steam.mesh.Mesh(copy=orig)

        xmin = cube.get_comp("1")
        xmax = cube.get_comp("2")
        ymin = cube.get_comp("3")
        ymax = cube.get_comp("4")
        zmin = cube.get_comp("5")
        zmax = cube.get_comp("6")

        # Mirror about x
        cube.transform(["m0x"])

        ### x sides should swap position
        check = np.all(cube.xyz_el.values[xmin,0] == orig.xyz_el.values[xmax,0])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[xmax,0] == orig.xyz_el.values[xmin,0])
        self.assertTrue( check )

        ### Other faces stay the same
        check = np.all(cube.xyz_el.values[ymin,1] == orig.xyz_el.values[ymin,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymax,1] == orig.xyz_el.values[ymax,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmin,2] == orig.xyz_el.values[zmin,2])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmax,2] == orig.xyz_el.values[zmax,2])
        self.assertTrue( check )


        # Mirror about non-zero y
        cube = steam.mesh.Mesh(copy=orig)
        yloc = 2.0
        cube.transform(["m{}y".format(yloc)])

        ### y sides should swap with offsets
        check = np.all(cube.xyz_el.values[ymin,1] == orig.xyz_el.values[ymin,1] + (yloc+0.5)*2.0)
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymax,1] == orig.xyz_el.values[ymax,1] + (yloc-0.5)*2.0)
        self.assertTrue( check )

        ### Other faces stay the same
        check = np.all(cube.xyz_el.values[xmin,0] == orig.xyz_el.values[xmin,0])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[xmax,0] == orig.xyz_el.values[xmax,0])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmin,2] == orig.xyz_el.values[zmin,2])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmax,2] == orig.xyz_el.values[zmax,2])
        self.assertTrue( check )

        # Mirror about non-zero y
        cube = steam.mesh.Mesh(copy=orig)
        zloc = -5.0
        cube.transform(["m{}z".format(zloc)])

        ### z sides should have offset
        check = np.all(cube.xyz_el.values[zmin,2] == orig.xyz_el.values[zmax,2] + 2*zloc)
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmax,2] == orig.xyz_el.values[zmin,2] + 2*zloc)
        self.assertTrue( check )

        ### Other faces rotate
        check = np.all(cube.xyz_el.values[ymin,1] == orig.xyz_el.values[ymin,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymax,1] == orig.xyz_el.values[ymax,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[xmin,0] == orig.xyz_el.values[xmin,0])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[xmax,0] == orig.xyz_el.values[xmax,0])
        self.assertTrue( check )

        
#! Run the tests
if __name__ == '__main__':
    ut.main()

