#! /bin/env python
# test_pointer_hashing.py

import unittest as ut
import steam

@ut.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh:")
class TestPointerHashing( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Create temporary files."""
        ### Create what's needed to make and store many solutions and meshes
        cls.n_mesh = 1000
        cls.mesh_list = []
        cls.lmesh_list = []
        cls.es_list = []

        ### Create a LibMeshInit object
        cls.init = steam.libmesh.Init()

        ### Make a bunch of meshes in python and libmesh
        for i in range(cls.n_mesh):
            mesh = steam.mesh.mk_cube( 1, 3, file = 0 )
            lmesh = steam.libmesh.Mesh.from_python( mesh )
            cls.mesh_list.append( mesh )
            cls.lmesh_list.append( lmesh )

    def test_mesh_storage( self ):
        """ Make a bunch of Meshes and store them in memory
        """
        ### Confirm that there are exactly the same number of meshes 
        ### and entries in the pointer storage dictionary
        ###     -- The "+ 1" in n_mesh + 1 is because of the Init object
        for lmesh in self.lmesh_list:
            self.assertTrue( len(lmesh.pointer_dict.keys()) == self.n_mesh + 1 )

            ### Make another mesh that points to the same libmesh object
            self.assertTrue( lmesh.pointer_dict[id(lmesh.p)] == 1 )
            throw_away = steam.libmesh.Mesh( pointer = lmesh.p )
            self.assertTrue( lmesh.pointer_dict[id(lmesh.p)] == 2 )
            throw_away2 = steam.libmesh.Mesh( pointer = lmesh.p )
            self.assertTrue( lmesh.pointer_dict[id(lmesh.p)] == 3 )
            del throw_away, throw_away2
            self.assertTrue( lmesh.pointer_dict[id(lmesh.p)] == 1 )

    def test_es_storage( self ):
        """ Make a bunch of EquationSystems and store them in memory
        """
        ### Make EquationSystems and confirm that there are now 2 * n_mesh
        ### entries in the pointer storage dictionary
        es_list = []
        for lmesh in self.lmesh_list:

            es = steam.libmesh.EquationSystems(lmesh)
            es_list.append( es )

        ### Again, the "+ 1" in n_mesh + 1 is because of the Init object
        self.assertTrue( len(es.pointer_dict.keys()) == 2 * self.n_mesh + 1 )

    @classmethod
    def tearDownClass( cls ):
        """Clean out temporary files."""
        cls.lmesh_list = None

#! Run the tests
if __name__ == '__main__':
    ut.main()
