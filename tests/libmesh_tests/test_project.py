#! /bin/env python

import unittest as ut
import steam

@ut.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestProject( ut.TestCase ):

    radius = 1.3

    @classmethod
    def setUpClass( cls ):
        """Create temporary files."""
        # Make grids and container
        cube = steam.mesh.mk_cube(2.0,21)
        sphere = steam.mesh.mk_sphere(cls.radius,31)
        cont = steam.container.Container("grids.h5")
        cont.add("/cube",cube)
        cont.add("/sphere",sphere)
        cont.write()

        ### Create a LibMeshInit object
        cls.init = steam.libmesh.Init()

    def test_project( self ):
        # Project the cube to the sphere
        cont = steam.container.Container("grids.h5")
        cont.read()
        cube   = steam.libmesh.Mesh.from_python(cont.obj['cube'])
        sphere = steam.libmesh.Mesh.from_python(cont.obj['sphere'])
        cube2sphere = cube.project(sphere)

        cube2 = steam.mesh.Mesh.from_libmesh(cube2sphere)
        # Check to make sure that all of the verticies have
        # the correct radius to within a percent
        for (i,j) in cube2.xyz_pt.iterrows():
            rad = (j['X']**2.0+j['Y']**2.0+j['Z']**2.0)**0.5
            #print(i,rad)
            error = abs(self.radius - rad)
            self.assertLess(error,self.radius*0.01)

    @classmethod
    def tearDownClass( cls ):
        """Clean out temporary files."""
        import os
        os.remove('grids.h5')
        #import shutil
        #shutil.rmtree('dir')


#! Run the tests
if __name__ == '__main__':
    ut.main()

