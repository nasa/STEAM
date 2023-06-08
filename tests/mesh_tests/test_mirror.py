#! /bin/env python
# test_mirror.py
""" Test the steam.mesh.half_to_full function
"""

import unittest as ut
import numpy as np
import steam
import random

from code import interact
#interact( local = dict( globals(), **locals() ) )

### Define planes for cutting along constant x-, y-, or z-planes
#   All planes should be right-handed
xplane = [(0,0,0), (0,1,0), (0,0,1)]        # x = 0
yplane = [(0,0,0), (0,0,1), (1,0,0)]        # y = 0
zplane = [(0,0,0), (1,0,0), (0,1,0)]        # z = 0

def coord_range( df ):
    """ Given a Pandas DataFrame, return a Series containing the range of
    each variable.
    """

    return df.max() - df.min()

class TestProject( ut.TestCase ):

    def standard_checks( self, full_mesh ):
        """ Simple suite of checks to run given a full body mesh made from
        self.cube.  This method is intended to be called by several or all of
        the other test methods
        """

        ### Full mesh should be a "Mesh"
        self.assertTrue( isinstance( full_mesh, steam.mesh.Mesh ) )

        ### The full mesh should have twice as many elements as the original
        self.assertEqual( full_mesh.conn.shape[0], 2*self.cube.conn.shape[0] )

        ### Output from coord_range should be double the half_mesh in one axis
        half_range = coord_range( self.cube.xyz_pt )
        full_range = coord_range( full_mesh.xyz_pt )
        self.assertEqual( len(full_range[full_range == half_range]), 2 )
        self.assertEqual( len(full_range[full_range != half_range]), 1 )
        self.assertEqual( len(full_range[full_range == 2 * half_range]), 1 )

    def setUp( self ):

        ### Create a half-body mesh from a segment of a cube
        self.cube_params = {'len': 5,
                            'npts' : 51}

        self.cube_params['x_offset'] =  random.uniform( 
                                        -self.cube_params['len'] /2.0,
                                         self.cube_params['len'] /2.0 )
        self.cube_params['y_offset'] =  random.uniform( 
                                        -self.cube_params['len'] /2.0,
                                         self.cube_params['len'] /2.0 )
        self.cube_params['z_offset'] =  random.uniform( 
                                        -self.cube_params['len'] /2.0,
                                         self.cube_params['len'] /2.0 )

        self.cube = steam.mesh.mk_cube( self.cube_params['len'], 
                                        self.cube_params['npts'] )

    def test_plane_input( self ):
        """ Confirm inputs and output type for mirror_plane argument. """

        good_inputs = ['x', 'y', 'z', 'X', 'Y', 'Z' ]
        bad_inputs = ['xx', 'xy', '1', 'a', 'B']
        list_inputs = [ ['X'], ['x', 'y']  ]

        ### Confirm that a ~steam.mesh.Mesh is returned
        for val in good_inputs:
            self.assertTrue( isinstance( steam.mesh.half_to_full( 
                                         self.cube, val ), steam.mesh.Mesh ) )
        for val in bad_inputs:
            with self.assertRaises( ValueError ):
                steam.mesh.half_to_full( self.cube, val )

        for val in list_inputs:
            with self.assertRaises( TypeError ):
                steam.mesh.half_to_full( self.cube, val )

    def test_mirror_0( self ):
        """ Shift the cube and mirror it about the zero planes """

        ### Shift the cube and cut it along the x, y, and z-planes, keeping
        ### all of the points in the positive octant
        self.cube.transform( [f't{self.cube_params["x_offset"]}X',
                              f't{self.cube_params["y_offset"]}Y',
                              f't{self.cube_params["z_offset"]}Z' ] )

        self.cube.slice( [xplane, yplane, zplane], 
                         delete=[[-1,1],[-1,1],[-1,1],] )
    
        ### Ensure we're starting in the right place -- there should be no 
        #   points with negative coordinates
        for coord in ('X', 'Y', 'Z'):
            self.assertAlmostEqual( self.cube.xyz_pt.min()[ coord ], 0 )

        for coord in ('X', 'Y', 'Z'):
            full_mesh = steam.mesh.half_to_full( self.cube, coord )
            self.standard_checks( full_mesh )
            self.assertEqual( full_mesh.xyz_pt.min()[ coord ], 
                                    -1.0 * self.cube.xyz_pt.max()[ coord ] )
        
    def test_on_plane_ind( self ):
        """ Test the steam.mesh._on_plane_indices function """

        ### By default the central axis of the sphere from mk_sphere is along
        #   the y-axis.  The z=0 plane always has points in it.  The x=0 plane
        #   will also always have points in it if the 'npts' input is even.

        for npts in range( 5, 101):
            sphere = steam.mesh.mk_sphere( 1, npts )

            on_x, off_x = steam.mesh._on_plane_indices( sphere, 'X', 0, 1e-12 )
            on_y, off_y = steam.mesh._on_plane_indices( sphere, 'Y', 0, 1e-12 )
            on_z, off_z = steam.mesh._on_plane_indices( sphere, 'Z', 0, 1e-12 )

            self.assertEqual( sphere.xyz_pt.shape[0], len(on_x) + len(off_x) )
            self.assertEqual( sphere.xyz_pt.shape[0], len(on_y) + len(off_y) )
            self.assertEqual( sphere.xyz_pt.shape[0], len(on_z) + len(off_z) )

            ### If npts is even, the number of points on the y=0 plane should 
            #   be 0.  Additionally, the number of points on the x=0 plane
            #   should be 2 -- the poles of the sphere.
            #   If npts is odd then the number of points on y=0  should be 
            #   equal to the number of points on the z=0 plane.
            if npts % 2 == 1:
                # Odd
                self.assertEqual( len(on_y), len(on_z) )
                self.assertEqual( len(off_y), len(off_z) )
            elif npts % 2 == 0:
                # Even
                self.assertEqual( len(on_x), 2 )
                self.assertEqual( len(on_y), 0 )

            ### Test the 'return_mirror_plane' input
            for coord, inds in zip( ('X', 'Y', 'Z'), (on_x, on_y, on_z) ):
                plane = steam.mesh.half_to_full( sphere, coord, mirror_coord=0.0, 
                                      return_mirror_plane=True )
                self.assertTrue( isinstance( plane, steam.mesh.Mesh ) )
                self.assertTrue( plane.is_point_cloud )
                self.assertTrue( plane.xyz_pt.equals( sphere.xyz_pt.loc[inds] ) )

#! Run the tests
if __name__ == '__main__':
    ut.main()

