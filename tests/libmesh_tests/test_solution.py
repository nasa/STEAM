#! /bin/env python

import unittest as ut
import steam
from filecmp import cmp
#from steam.libmesh.equation_systems import ElemNodeStorage as DataAt

@ut.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh:")
class TestSolution( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Create temporary files."""
        ### Set up solution and mesh
        #   Uniform cube
        cls.mesh = steam.mesh.mk_cube( 1, 3, file = 2 )
        cls.soln = steam.solution.uniform_soln( cls.mesh, 
                    node_elem_flag = 'ELEMENTS', values = ['rand']*20 )

        ### Create a LibMeshInit object
        cls.init = steam.libmesh.Init()

    def test_triq_to_libmesh( self ):
        """Test triq -> libmesh."""

        lmesh = steam.libmesh.Mesh.from_python(self.mesh)

        es = steam.libmesh.EquationSystems(lmesh)
        es.from_python(self.soln)
        es.write('PLT','test.dat')
        es.write('E','test.e')

    def test_soln_from_libmesh_ele( self ):
        """Get the solution back from libMesh.
        Solution at elements.

        NOTE: This test is reproduced by ../solution_tests/test_conversion.py
        """

        ### Send the solution to libMesh
        lmesh = steam.libmesh.Mesh.from_python(self.mesh)
        es = steam.libmesh.EquationSystems(lmesh)
        sys_name = es.from_python(self.soln, soln_name = 'test_sys' )

        ### Get it back
        new_soln = steam.solution.Solution.from_libmesh( es, name = sys_name )

        ### Compare the data
        self.assertTrue( new_soln.data.equals(self.soln.data), 
                         "new_soln.data does not equal soln.data" )
        
        ### Write the two solutions as triq files and compare them
        original_out_file = 'original_soln.triq'
        new_out_file      = 'new_soln.triq'
        steam.io.write_triq_ascii( soln = self.soln, name = original_out_file)
        steam.io.write_triq_ascii( soln = new_soln, name = new_out_file)
        self.assertTrue( cmp(original_out_file, new_out_file, shallow = False))

    def test_soln_from_libmesh_nodes( self ):
        """Get the solution back from libMesh.
        Solution at nodes.

        NOTE: This test is reproduced by ../solution_tests/test_conversion.py
        """

        ### Send the solution to libMesh
        lmesh = steam.libmesh.Mesh.from_python(self.mesh)
        es = steam.libmesh.EquationSystems(lmesh)
        node_soln = self.soln.element_to_node( return_new = True )
        sys_name = es.from_python(node_soln, soln_name = 'nodes' )

        ### Get it back
        new_soln = steam.solution.Solution.from_libmesh( es, name = sys_name )

        ### Compare the data
        self.assertTrue( new_soln.data.equals(node_soln.data), 
                         "new_soln.data does not equal node_soln.data" )

        ### Write the two solutions as triq files and compare them
        original_out_file = 'original_soln.triq'
        new_out_file      = 'new_soln.triq'
        steam.io.write_triq_ascii( soln = node_soln, name = original_out_file)
        steam.io.write_triq_ascii( soln = new_soln, name = new_out_file)
        self.assertTrue( cmp(original_out_file, new_out_file, shallow = False))

    @classmethod
    def tearDownClass( cls ):
        """Clean out temporary files."""
        import os
        files_to_remove = ( 'test.dat', 'test.e', 'original_soln.triq', 
                            'new_soln.triq', 'cube.tri' )
        for f in files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )

#! Run the tests
if __name__ == '__main__':
    ut.main()
