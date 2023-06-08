#!/usr/bin/env python

import unittest
import steam

@unittest.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestInit( unittest.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """ This method runs first and exactly once, before all other methods.
        """

        ### Initialize libmesh
        cls.init = steam.libmesh.Init()
        pass

    def test_delete(self):
        """ This test is going to delete out of order.

        This should ensure that we're robust with regards to garbage collection.
        """

        cube = steam.mesh.mk_cube(5.0,31)
        mesh = steam.libmesh.Mesh.from_python(cube)
        es   = steam.libmesh.EquationSystems(mesh)
        es.add_transient_explicit_system("MyTransientExplicitSystem")
        del mesh
        del es

### Run the tests
if __name__ == '__main__':
    unittest.main()
