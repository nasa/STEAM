#! /bin/env python

import math
import numpy as np
import os
import pandas as pd
import unittest as ut

import steam

class TestProject( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Make a simple cube test case to work with.
        """

        ### Save cube parameters for use in other tests
        cls.cube = {}
        cls.cube['length'] = 1
        cls.cube['dims']   = 5

        mesh = steam.mesh.mk_cube( cls.cube['length'], cls.cube['dims'])
        mesh.set_comp("One"    ,["1"]        )
        mesh.set_comp("Even"   ,["2","4","6"]    )
        mesh.set_comp("Odd"    ,["One","3","5"])
        mesh.set_comp("First5" ,[1,2,3,4,5])

        cls.mesh_file = 'test.h5'
        with steam.container.Container( cls.mesh_file, 'w' ) as cont:
            cont.add("/cube",mesh)

        cls.files_to_remove = []
        cls.files_to_remove.append( cls.mesh_file )

    def test_point_cloud( self ):
        """Make sure that a mesh is converted to a point cloud properly
        """

        cont = steam.container.Container("test.h5")
        cont.read()
        cube = cont.obj['cube']

        self.assertFalse( cube.is_point_cloud )
        self.assertTrue( cube.static_soln is None )

        ### Convert to point cloud
        cube.to_point_cloud()

        self.assertTrue( cube.is_point_cloud )
        self.assertTrue( cube.xyz_pt.equals( cube.xyz_el ) )

    def test_read( self ):
        """Make sure reads are lossless.
        """
        mesh = steam.mesh.mk_cube( 1, 5)
        mesh.set_comp("One" ,["1"]        )
        mesh.set_comp("Even",["2","4","6"]    )
        mesh.set_comp("Odd" ,["One","3","5"])

        cont = steam.container.Container("test.h5")
        cont.read()
        cube = cont.obj['cube']

        #! Make sure that each component is identical to what was just made
        for comp in mesh.comp_table:
            new_list = mesh.get_comp(comp)
            old_list = cube.get_comp(comp)

            self.assertEqual(new_list,old_list)

    def test_elem( self ):
        """Check to make sure components of elements work.
        """

        cont = steam.container.Container("test.h5")
        cont.read()
        cube = cont.obj['cube']

        comp = ["First5"]
        new_list = [1,2,3,4,5]
        old_list = cube.get_comp(comp)

        self.assertEqual(new_list,old_list)

    def test_numbers( self ):
        """Check to make sure numerical components work.
        """

        cont = steam.container.Container("test.h5")
        cont.read()
        cube = cont.obj['cube']

        #! Make sure that each component is identical to what was just made
        for comp in [1,2,3,4,5,6]:
            subset   = cube.comp[cube.comp[0] == comp]
            new_list = list(subset.index)
            old_list = cube.get_comp([comp])

            self.assertEqual(new_list,old_list)

    def test_strings( self ):
        """Check to make sure string components work.
        """

        cont = steam.container.Container("test.h5")
        cont.read()
        cube = cont.obj['cube']

        comp = "One"
        old_list = cube.get_comp([comp])
        new_list = list(cube.comp[cube.comp[0] == 1].index)
        self.assertEqual(new_list,old_list)

        comp = "Odd"
        old_list = cube.get_comp([comp])
        new_list.extend(list(cube.comp[cube.comp[0] == 3].index))
        new_list.extend(list(cube.comp[cube.comp[0] == 5].index))
        self.assertEqual(new_list,old_list)

    def test_comp_by_xyz( self ):
        """Check to make sure the 'set_comp_by_xyz' works.
        """

        ### Create the mesh
        radius = 10
        npts   = 11
        mesh = steam.mesh.mk_sphere( radius, npts )
        mesh.get_xyz_el()

        ### Define the component ranges
        xlim = (6.0, 8.0)        # (xmin, xmax)
        ylim = (0.0, radius*2)   # (ymin, ymax)
        zlim = (0.0, radius*2)   # (zmin, zmax)

        ### We'll use DataFrame.query to determine the correct elements in comp
        query_string = ( f'{xlim[0]} < X and X < {xlim[1]} and '
                       + f'{ylim[0]} < Y and Y < {ylim[1]} and '
                       + f'{zlim[0]} < Z and Z < {zlim[1]} ' )
        correct_elems = mesh.xyz_el.query( query_string ).index.to_list()

        mesh.set_comp_by_xyz('ring', xmin=xlim[0], xmax=xlim[1], 
                                     ymin=ylim[0], ymax=ylim[1], 
                                     zmin=zlim[0], zmax=zlim[1] )
        ### Element order may vary, so compare sets instead of lists
        self.assertEqual( set(mesh.get_comp(['ring'])), set(correct_elems) )

    def test_comp_by_x_and_theta( self ):
        """Check to make sure the 'set_comp_by_x_and_theta' works.
        """

        ### NOTE: Mesh.set_comp_by_x_and_theta is very poorly documented and 
        #   seems to fill a niche role.  As a result, Cam Ebner edited it to
        #   remove all functionality on 3/31/21.  The only thing it should
        #   ever do is raise a NotImplementedError when called.

        ### Create the mesh
        radius = 10
        npts   = 11
        mesh = steam.mesh.mk_sphere( radius, npts )

        with self.assertRaises( NotImplementedError ):
            mesh.set_comp_by_x_and_theta( 'ring' )
#        mesh.get_xyz_el()
#
##        mesh = steam.mesh.mk_sphere (10.0, 11)
#        # tmax is greater than tmin
#        correct_elem_numbers = [0, 1, 2, 18, 19, 20, 36, 37, 38, 54, 55, 56, 73, 74, 307, 308, 324, 325, 326, 342, 343, 344, 360, 361, 362]
#        mesh.set_comp_by_x_and_theta ('ring', xmin=1.0, xmax=5.0, tmin=60.0, tmax=120.0)
#        self.assertEqual (mesh.comp_table ['ring']['elem'], correct_elem_numbers)
#        # tmin is greater than tmax (bridges theta = 0 line)
#        correct_elem_numbers = [47, 61, 62, 63, 64, 66]
#        mesh.set_comp_by_x_and_theta ('ring2', xmin=6.0, xmax=8.0, tmin=300.0, tmax=30.0)
#        self.assertEqual (mesh.comp_table ['ring2']['elem'], correct_elem_numbers)

    def test_add_comps_after_read( self ):
        """ Read the pre-made mesh and add components to it.
        """
        with steam.container.Container( self.mesh_file, 'r' ) as cont:
            cont.read()
        mesh = cont.obj['cube']

        ### Confirm that we can write the mesh to a new h5 file
        new_out = 'same_mesh.h5'
        self.files_to_remove.append( new_out )
        with steam.container.Container( new_out, 'w' ) as cont2:
            cont2.add( 'mesh', mesh )
        self.assertTrue( os.path.isfile( new_out ) )

        ### Load the mesh from the new output
        with steam.container.Container( new_out, 'r' ) as in_cont:
            in_cont.read()
        self.assertEqual( mesh, in_cont.obj['mesh'] )
        mesh = in_cont.obj['mesh']

        ### Add another component and confirm that we can still write the mesh
        mesh.set_comp("Second5", [6,7,8,9,10] )
        second_out = 'one_new_comp.h5'
        self.files_to_remove.append( second_out )
        with steam.container.Container( second_out, 'w' ) as cont3:
            cont3.add( 'mesh', mesh )
        self.assertTrue( os.path.isfile( second_out ) )

        ### Add a third component and confirm that we can still write the mesh
        mesh.set_comp("new_comps", ['Second5', 'Third5'] )
        third_out = 'two_new_comp.h5'
        self.files_to_remove.append( third_out )
        with steam.container.Container( third_out, 'w' ) as cont4:
            cont4.add( 'mesh', mesh )
        self.assertTrue( os.path.isfile( third_out ) )

    def test_huge_component( self ):
        """ Create a component with a million entries and try to save the mesh
        """

        with steam.container.Container( self.mesh_file, 'r' ) as cont:
            cont.read()
        mesh = cont.obj['cube']

        mesh.set_comp("gigantic", np.ones(int(1e6)) )

        ### Attempt to write the mesh with large component
        out_cont = 'big_comps.h5'
        self.files_to_remove.append( out_cont )
        with steam.container.Container( out_cont, 'w' ) as out:
            out.add( 'mesh', mesh )
        self.assertTrue( os.path.isfile( out_cont ) )

    def test_old_comp_table( self ):
        """ Ensure backwards compatibility with old comp_table location
        """

        old_cont = 'old_comp_table.h5'

        ### Manually open the hdf5 file to ensure the comp_table is stored
        ### in the old fashion
        #   This line copied and altered from container.open_store:
        hdf5 = pd.HDFStore( old_cont, complevel=9, complib='blosc', mode='r' )
        #   "iroot" must exist for reading hdf5
        root = f"/{steam.util.string_to_hdf5_key( 'mesh' )}"
        iroot = f'{root}/xyz_pt/table'
        with self.assertRaises(TypeError):
            comp_table = steam.util.read_hdf5_pnode( hdf5,
                                                     f'{root}/comp_table' )
        ### Newer comp_table storage location should fail
        comp_table = hdf5.get_storer(iroot).attrs.comps
        #   This line copied and altered from container.close_store:
        hdf5.close()

        with steam.container.Container( old_cont, 'r' ) as cont:
            cont.read()

        old_mesh = cont.obj['mesh']

        ### Read newly created mesh to compare with the old one
        with steam.container.Container( self.mesh_file, 'r' ) as new_cont:
            new_cont.read()
            new_mesh = new_cont.obj['cube']

        self.assertEqual( old_mesh, new_mesh )
        self.assertEqual( new_mesh.comp_table, comp_table )
        self.assertEqual( old_mesh.comp_table, comp_table )

    def test_simple_to_mesh( self ):
        """ Test the component allocation in mesh.simple_to_mesh
        """

        ### Get the mesh created in setUpClass
        with steam.container.Container( self.mesh_file, 'r' ) as cont:
            cont.read()
        original_cube = cont.obj['cube']

        ### Construct the same cube semi-manually
        nodes, elems, comps = steam.mesh.build_cube( self.cube['length'], 
                                                     self.cube['dims'] )
        same_cube = steam.mesh.simple_to_mesh( nodes, elems, comps )

        ### same_cube is not exactly equal to original_cube because
        ### original_cube had string components added
        self.assertTrue( original_cube.xyz_pt.equals(same_cube.xyz_pt) )
        self.assertTrue( original_cube.conn.equals(same_cube.conn) )
        original_cube.get_xyz_el()
        same_cube.get_xyz_el()
        self.assertTrue( original_cube.xyz_el.equals(same_cube.xyz_el) )

        ### simple_to_mesh should fail if given a component that's a string
        string_comp = list( comps )
        string_comp[-1] = 'string'
        with self.assertRaises( TypeError ):
            string_cube = steam.mesh.simple_to_mesh(nodes, elems, string_comp)

        ### simple_to_mesh should fail if given a component list 
        ### that's too long
        long_comp = list( comps )
        long_comp.append( 1 )
        with self.assertRaises( ValueError ):
            string_cube = steam.mesh.simple_to_mesh(nodes, elems, long_comp)

    @classmethod
    def tearDownClass( cls ):
        """Clean up.
        """

        for f in cls.files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )


#! Run the tests
if __name__ == '__main__':
    ut.main()

