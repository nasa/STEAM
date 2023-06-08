 #! /bin/env python

import unittest as ut
import numpy as np
import steam

class TestProject( ut.TestCase ):



    def test_translation( self ):
        """ """
        ### Make the cubes
        orig = steam.mesh.mk_cube( 1, 5)
        orig.get_xyz_el()
        cube = steam.mesh.Mesh(copy=orig)

        xdist = -1.0
        ydist =  2.0
        zdist = 10.0

        # One at a time
        cube.transform(["t{}x".format(xdist)])

        ### Check
        check = np.all(cube.xyz_pt.values[:,0] - orig.xyz_pt.values[:,0] == xdist)
        self.assertTrue( check )
        
        # Two at a time
        cube.transform(["t{}y".format(ydist),"t{}z".format(zdist)])

        ### Check
        check = np.all(cube.xyz_pt.values[:,1] - orig.xyz_pt.values[:,1] == ydist)
        self.assertTrue( check )

        check = np.all(cube.xyz_pt.values[:,2] - orig.xyz_pt.values[:,2] == zdist)
        self.assertTrue( check )
        
    def test_90deg_rotation( self ):
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

        # Do x rotation
        cube.transform(["r90x"])

        ### Check sides 1/2 (which should not change with x rotation)
        check = np.all(cube.xyz_el.values[xmin,0] == orig.xyz_el.values[xmin,0])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[xmax,0] == orig.xyz_el.values[xmax,0])
        self.assertTrue( check )

        ### Other faces rotate
        check = np.all(cube.xyz_el.values[zmin,1] == orig.xyz_el.values[ymax,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmax,1] == orig.xyz_el.values[ymin,1])
        self.assertTrue( check )

        check = np.all(cube.xyz_el.values[ymin,2] == orig.xyz_el.values[zmin,2])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymax,2] == orig.xyz_el.values[zmax,2])
        self.assertTrue( check )


        # Do y rotation
        cube = steam.mesh.Mesh(copy=orig)
        cube.transform(["r90y"])

        ### Check sides 1/2 (which should not change with x rotation)
        check = np.all(cube.xyz_el.values[ymin,1] == orig.xyz_el.values[ymin,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymax,1] == orig.xyz_el.values[ymax,1])
        self.assertTrue( check )

        ### Other faces rotate
        check = np.all(cube.xyz_el.values[zmin,0] == orig.xyz_el.values[xmin,0])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmax,0] == orig.xyz_el.values[xmax,0])
        self.assertTrue( check )

        check = np.all(cube.xyz_el.values[xmin,2] == orig.xyz_el.values[zmax,2])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[xmax,2] == orig.xyz_el.values[zmin,2])
        self.assertTrue( check )

        # Do z rotation
        cube = steam.mesh.Mesh(copy=orig)
        cube.transform(["r90z"])

        ### Check sides 1/2 (which should not change with x rotation)
        check = np.all(cube.xyz_el.values[zmin,2] == orig.xyz_el.values[zmin,2])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[zmax,2] == orig.xyz_el.values[zmax,2])
        self.assertTrue( check )

        ### Other faces rotate
        check = np.all(cube.xyz_el.values[xmin,1] == orig.xyz_el.values[ymin,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[xmax,1] == orig.xyz_el.values[ymax,1])
        self.assertTrue( check )

        check = np.all(cube.xyz_el.values[ymin,0] == orig.xyz_el.values[xmax,0])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymax,0] == orig.xyz_el.values[xmin,0])
        self.assertTrue( check )


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

        
##
    def test_comp_transform( self ):
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
        
        # translate ymin in x
        xdist = 5.0
        comps = ["3"]
        cube.transform(["t{}x".format(xdist)],comps)

        ### Check side 3 (which should change with x translation)
        check = np.all(cube.xyz_el.values[ymin,0] - orig.xyz_el.values[ymin,0] == xdist)
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymin,1] == orig.xyz_el.values[ymin,1])
        self.assertTrue( check )
        check = np.all(cube.xyz_el.values[ymin,2] == orig.xyz_el.values[ymin,2])
        self.assertTrue( check )

        ### Check side 4 (which should not have changed with x translation)
        check = np.all(cube.xyz_el.values[ymax,:] == orig.xyz_el.values[ymax,:])
        self.assertTrue( check )
        
        ### Make the cube
        cube= steam.mesh.Mesh(copy=orig)
        
        xmin = cube.get_comp("1")
        xmax = cube.get_comp("2")
        ymin = cube.get_comp("3")
        ymax = cube.get_comp("4")
        zmin = cube.get_comp("5")
        zmax = cube.get_comp("6")
        
        # translate ymin in x
        xdist = 2.5
        comps = ["2","4"]
        cube.transform(["t{}x".format(xdist),"r180x"],comps)
        #cube_2.transform(["r180x"],comps_2)
        
        ### Check side 2 and 4 (which should change with x translation and rotation about x)
        
        
        
        check = np.all(cube.xyz_el.values[xmax,0] == orig.xyz_el.values[xmin,0] + xdist + 1.0)
        self.assertTrue( check )
        #check = np.all(cube_2.xyz_el.values[xmax,1] + orig.xyz_el.values[xmin,1] > 1.0)
        #self.assertTrue( check )
        #check = np.all(cube_2.xyz_el.values[xmax,2] == orig.xyz_el.values[xmin,2])
        #self.assertTrue( check )
        
        check = np.all(cube.xyz_el.values[ymax,0] == orig.xyz_el.values[ymin,0] + xdist)
        self.assertTrue( check )
        #check = np.all(cube_2.xyz_el.values[ymax,1] == orig.xyz_el.values[ymin,1])
        #self.assertTrue( check )
        #check = np.all(cube_2.xyz_el.values[ymax,2] == orig.xyz_el.values[ymin,2])
        #self.assertTrue( check )
        
#! Run the tests
if __name__ == '__main__':
    ut.main()

