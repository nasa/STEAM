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
        cls.prime_cube =    steam.mesh.mk_cube( 1, cube_pts)
        cls.prime_sphere =  steam.mesh.mk_sphere( 1, sphere_pts)

        ### Set up a static solution for the sphere
        cls.prime_sphere.static_soln = steam.solution.uniform_soln(
                           cls.prime_sphere, 'ELEM', [1.25, 'rand', 3.141] )

    def test_copy( self ):
        """Copy each mesh and confirm that the copies equal originals"""

        copy_cube = steam.mesh.Mesh( copy = self.prime_cube )
        copy_sphere = steam.mesh.Mesh( copy = self.prime_sphere )

        ### Cube
        self.assertFalse( copy_cube is self.prime_cube )
        self.assertEqual( self.prime_cube, copy_cube )

        ### Sphere
        self.assertFalse( copy_sphere is self.prime_sphere )
        self.assertEqual( self.prime_sphere, copy_sphere )

    def test_new( self ):
        """Make new meshes that are the same as the prime meshes and compare"""

        new_cube =    steam.mesh.mk_cube( 1, cube_pts)
        new_sphere =  steam.mesh.mk_sphere( 1, sphere_pts)

        ### Cube
        self.assertFalse( new_cube is self.prime_cube )
        self.assertEqual( self.prime_cube, new_cube )

        ### Sphere
        self.assertFalse( new_sphere is self.prime_sphere )
        self.assertTrue( steam.mesh.Mesh.__eq__( self.prime_sphere, new_sphere,
                         compare_static = False ) )

    def test_flags( self ):
        """Test the functionality of the flags"""

        new_cube =    steam.mesh.mk_cube( 1, cube_pts)
        new_sphere =  steam.mesh.mk_sphere( 1, sphere_pts)

        ### Test compare_static -- can only be dones with the sphere
        self.assertNotEqual( self.prime_sphere, new_sphere )
        #   Use __eq__
        self.assertFalse( steam.mesh.Mesh.__eq__(self.prime_sphere, new_sphere,
            compare_static=True) )
        self.assertTrue( steam.mesh.Mesh.__eq__(self.prime_sphere, new_sphere,
            compare_static=False) )
        #   Use __neq__
        self.assertTrue( steam.mesh.Mesh.__neq__(self.prime_sphere, new_sphere,
            compare_static=True) )
        self.assertFalse( steam.mesh.Mesh.__neq__(self.prime_sphere, new_sphere,
            compare_static=False) )
        #   Copy the static solution and compare
        new_sphere.static_soln = steam.solution.Solution( 
                                    copy = self.prime_sphere.static_soln )
        self.assertEqual( new_sphere, self.prime_sphere )
        self.assertTrue( steam.mesh.Mesh.__eq__(self.prime_sphere, new_sphere,
            compare_static=True) )
        self.assertFalse( steam.mesh.Mesh.__neq__(self.prime_sphere, new_sphere,
            compare_static=True) )

        ### Test compare_xyz_el on the cube
        self.assertEqual( self.prime_cube, new_cube )
        new_cube.get_xyz_el()
        self.assertEqual( self.prime_cube, new_cube )
        #   Use __eq__
        self.assertFalse( steam.mesh.Mesh.__eq__(self.prime_cube, new_cube,
            compare_xyz_el=True) )
        self.assertTrue( steam.mesh.Mesh.__eq__(self.prime_cube, new_cube,
            compare_xyz_el=False) )
        #   Use __neq__
        self.assertTrue( steam.mesh.Mesh.__neq__(self.prime_cube, new_cube,
            compare_xyz_el=True) )
        self.assertFalse( steam.mesh.Mesh.__neq__(self.prime_cube, new_cube,
            compare_xyz_el=False) )

    def test_hash( self ):
        """Test the mesh.update_hash function"""

        identical_cube = steam.mesh.mk_cube( 1, cube_pts)
        self.assertEqual( self.prime_cube, identical_cube )

        ### Hashes may or may not be populated
        self.prime_cube.update_hash()
        identical_cube.update_hash()
        self.assertEqual( self.prime_cube.hash, identical_cube.hash )

        ### Meshes should no longer be equal once one has a static_soln
        identical_cube.static_soln = steam.solution.uniform_soln(
            identical_cube, values=['rand', 'rand'] )
        identical_cube.update_hash()
        self.assertNotEqual( self.prime_cube.hash, identical_cube.hash )
        self.assertNotEqual( self.prime_cube, identical_cube )

#! Run the tests
if __name__ == '__main__':
    ut.main()

