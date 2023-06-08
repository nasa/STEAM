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

#        self.cube_params['x_offset'] =  2.2538706455661934
#        self.cube_params['y_offset'] =  1.103165596678016
#        self.cube_params['z_offset'] =  2.442885688382634
        self.cube = steam.mesh.mk_cube( self.cube_params['len'], 
                                        self.cube_params['npts'] )

        self.centerline_tol = 1e-8

#    @ut.skip( 'Skip for running the whole test suite.' )
    def test_mirror_0( self ):
        """ Shift the cube and mirror it about the zero planes """

        ### Shift the cube and cut it along the x, y, and z-planes, keeping
        ### all of the points in the positive octant
        self.cube.transform( [f't{self.cube_params["x_offset"]}X',
                              f't{self.cube_params["y_offset"]}Y',
                              f't{self.cube_params["z_offset"]}Z' ] )

        self.cube.slice( [xplane, yplane, zplane], 
                         delete=[[-1,1],[-1,1],[-1,1],] )

        ### Create the solution for what remains of the cube
        self.half_elem_soln = steam.solution.uniform_soln( self.cube, 'ELEM',
                                            values=['rand','rand'] )
        self.half_node_soln = steam.solution.uniform_soln( self.cube, 'NODE',
                                            values=['rand','rand'] )
        self.assertTrue( 
                self.cube.xyz_el.index.equals( self.half_elem_soln.data.index ) )
    
        for coord in ('X', 'Y', 'Z'):
            ### Ensure we're starting in the right place -- there should be no 
            #   points with negative coordinates
            self.assertAlmostEqual( self.cube.xyz_pt.min()[ coord ], 0 )

            ### Add coordinate data to the solution for checking later
            self.half_elem_soln.add_var( coord, self.cube.xyz_el[coord], 
                                         index_in=self.cube.xyz_el.index )
            self.half_node_soln.add_var( coord, self.cube.xyz_pt[coord], 
                                         index_in=self.cube.xyz_pt.index )

            ### Create the full body mesh and solutions
            full_mesh = steam.mesh.half_to_full( self.cube, coord, 
                                                 tol=self.centerline_tol )
            elem_soln = steam.solution.half_to_full( full_mesh, 
                                                     self.half_elem_soln )
            node_soln = steam.solution.half_to_full( full_mesh, 
                                                     self.half_node_soln )

            ### Confirm that the mirrored data matches the unmirrored data
            #   Element solution
            self.assertTrue( np.array_equal( 
                    elem_soln.data.loc[elem_soln.mesh._elem_map.index].values, 
                    elem_soln.data.loc[elem_soln.mesh._elem_map.values].values 
                    ) )
            #   Node solution
            self.assertTrue( np.array_equal( 
                    node_soln.data.loc[node_soln.mesh._node_map.index].values, 
                    node_soln.data.loc[node_soln.mesh._node_map.values].values 
                    ) )

            ### Check that all coordinate data stored in solution is <=0 and 
            #   accurate
            full_mesh.get_xyz_el()
#            interact( local = dict( globals(), **locals() ) )
            self.assertTrue( np.all( elem_soln.data[coord] >= 
                                                    -1*self.centerline_tol ) )
            self.assertTrue( np.all( node_soln.data[coord] >= 
                                                    -1*self.centerline_tol ) )
            self.assertFalse( np.all( full_mesh.xyz_pt[coord] >= 0 ) )
            self.assertFalse( np.all( full_mesh.xyz_el[coord] >= 0 ) )

            #   Same side as half-mesh
            self.assertTrue( 
                elem_soln.data[coord].loc[full_mesh._elem_map.index].equals(
                full_mesh.xyz_el[coord].loc[full_mesh._elem_map.index] ) )
            self.assertTrue( 
                node_soln.data[coord].loc[full_mesh._node_map.index].equals(
                full_mesh.xyz_pt[coord].loc[full_mesh._node_map.index] ) )
            #   Mirrored side -- Note use of allclose because of minor 
            #   differences after full_mesh.get_xyz_el()
            self.assertTrue( np.allclose(
                elem_soln.data[coord].loc[full_mesh._elem_map.values],
                -1.0 * full_mesh.xyz_el[coord].loc[full_mesh._elem_map.values],
                atol=1e-12) )
            self.assertTrue( np.allclose(
                node_soln.data[coord].loc[full_mesh._node_map.values],
                -1.0 * full_mesh.xyz_pt[coord].loc[full_mesh._node_map.values],
                atol=1e-12) )

    def test_flip_vars( self ):
        """Exercise functionality of the flip_vars input """

        ### Shift the cube and cut it along the x, y, and z-planes, keeping
        ### all of the points in the positive octant
        self.cube.transform( [f't{self.cube_params["x_offset"]}X',
                              f't{self.cube_params["y_offset"]}Y',
                              f't{self.cube_params["z_offset"]}Z' ] )

        self.cube.slice( [xplane, yplane, zplane], 
                         delete=[[-1,1],[-1,1],[-1,1],] )

        ### Create the solution for what remains of the cube
        self.half_elem_soln = steam.solution.uniform_soln( self.cube, 'ELEM',
                                            values=['rand','rand'] )
        self.half_node_soln = steam.solution.uniform_soln( self.cube, 'NODE',
                                            values=['rand','rand'] )
        self.assertTrue( 
                self.cube.xyz_el.index.equals( self.half_elem_soln.data.index ) )
    
        for coord in ('X', 'Y', 'Z'):
            ### Ensure we're starting in the right place -- there should be no 
            #   points with negative coordinates
            self.assertAlmostEqual( self.cube.xyz_pt.min()[ coord ], 0 )

            ### Add coordinate data to the solution for checking later
            self.half_elem_soln.add_var( coord, self.cube.xyz_el[coord], 
                                         index_in=self.cube.xyz_el.index )
            self.half_node_soln.add_var( coord, self.cube.xyz_pt[coord], 
                                         index_in=self.cube.xyz_pt.index )

            ### Create the full body mesh and solutions
            full_mesh = steam.mesh.half_to_full( self.cube, coord,
                                                 tol=self.centerline_tol )
            full_mesh.get_xyz_el()
            elem_soln = steam.solution.half_to_full( full_mesh, 
                                         self.half_elem_soln, [coord] )
            node_soln = steam.solution.half_to_full( full_mesh, 
                                         self.half_node_soln, [coord] )
    
            self.assertTrue( node_soln.data[coord].equals( 
                             full_mesh.xyz_pt[coord] ) )
            ### Small deviations may occur in the data, so strict equality
            #   will not be enforced
            threshold = 1e-16
            self.assertTrue( (elem_soln.data[coord] - 
                              full_mesh.xyz_el[coord]).max() < threshold )

#! Run the tests
if __name__ == '__main__':
    ut.main()

