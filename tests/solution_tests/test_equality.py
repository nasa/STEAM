#! /bin/env python
# test_static_soln.py

import unittest as ut
import steam
#from code import interact
#interact( local = dict( globals(), **locals() ) )

cube_pts = 15
sphere_pts = 36

class TestEqualilty( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Make a simple cube test case to work with.
        """

        ### Cubes with solution at NODES
        #   Define a mesh that DOES NOT have xyz_el populated
        cls.mesh = steam.mesh.mk_cube( 1, cube_pts)
        #cls.prime_sphere =  steam.mesh.mk_sphere( 1, sphere_pts)

        ### Set up a static solution for the sphere
        cls.elem_soln = steam.solution.uniform_soln(
                           cls.mesh, 'ELEM', [1.25, 3.141] )
        cls.node_soln = steam.solution.uniform_soln(
                           cls.mesh, 'NODE', [2.01, 6.282] )

    def test_copy( self ):
        """Copy each mesh and confirm that the copies equal originals"""

        copy_elem = steam.solution.Solution( copy = self.elem_soln )
        copy_node = steam.solution.Solution( copy = self.node_soln )

        ### Element-based solution
        self.assertFalse( copy_elem is self.elem_soln )
        self.assertEqual( self.elem_soln, copy_elem )
        self.assertTrue( self.elem_soln == copy_elem )

        ### Node-based solution
        self.assertFalse( copy_node is self.node_soln )
        self.assertEqual( self.node_soln, copy_node )
        self.assertTrue( self.node_soln == copy_node )

    def test_new( self ):
        """Make new solutions that are the same as the original solutions 
        and compare.
        """

        new_mesh = steam.mesh.mk_cube( 1, cube_pts)
        new_elem = steam.solution.uniform_soln(new_mesh, 'ELEM', [1.25, 3.141])
        new_node = steam.solution.uniform_soln(new_mesh, 'NODE', [2.01, 6.282])

        ### Element-based solution
        self.assertFalse( new_elem is self.elem_soln )
        self.assertEqual( self.elem_soln, new_elem )
        self.assertTrue( self.elem_soln == new_elem )

        ### Node-based solution
        self.assertFalse( new_node is self.node_soln )
        self.assertEqual( self.node_soln, new_node )
        self.assertTrue( self.node_soln == new_node )

    def test_neq( self ):
        """Make new solutions that are the close to the original solutions
        and confirm that they don't qualify as equal.
        """

        new_node = steam.solution.uniform_soln( self.mesh, 'NODE', 
                                                [2.01, 6.282])

        ### Element-based solution
        self.assertNotEqual( self.elem_soln, 
                steam.solution.uniform_soln(self.mesh, 'ELEM', [1.25, 3.14]) )
        self.assertNotEqual( self.elem_soln, 
                steam.solution.uniform_soln(self.mesh, 'NODE', [1.25, 3.141]) )
        self.assertTrue( self.elem_soln !=
                steam.solution.uniform_soln(self.mesh, 'NODE', [1.25, 3.141]) )

        ### Other checks
        self.assertNotEqual( self.elem_soln, 'string' )
        empty_node = steam.solution.Solution()
        empty_node.store_at = 'NODE'
        empty_elem = steam.solution.Solution()
        empty_elem.store_at = 'ELEM'
        self.assertNotEqual( empty_node, empty_elem )
        self.assertTrue( empty_node != empty_elem )
        self.assertNotEqual( empty_node, self.node_soln )
        self.assertNotEqual( self.elem_soln, empty_elem )

#! Run the tests
if __name__ == '__main__':
    ut.main()

