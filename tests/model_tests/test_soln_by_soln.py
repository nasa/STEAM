#! /bin/env python

""" Replicate the LibMESH tests for the Python-side modeling methods.
"""

import unittest as ut
import steam
import os
from code import interact

#@ut.skip( 'Not implemented yet' )
class TestModels( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        #self.mesh.write_tri_ascii( 'cube' )
        cls.remove_files = []

    def setUp( self ):
        """Create a mesh and several solutions with which to test
        """

        ### Create a cube and populate it with uniform solutions
        self.mesh = steam.mesh.mk_cube( 1, 5, file = 0 )

        ### Two node-based and two elem-based solutions are needed for 
        #   soln-by-soln modeling
        self.node_solnA = steam.solution.uniform_soln( self.mesh, 
                          node_elem_flag = 'NODE', values = [1.0,'rand'] )
        self.node_solnB = steam.solution.uniform_soln( self.mesh, 
                          node_elem_flag = 'NODE', values = [2.0,'rand'] )
        self.elem_solnA = steam.solution.uniform_soln( self.mesh, 
                          node_elem_flag = 'ELEM', values = [1.0,'rand'] )
        self.elem_solnB = steam.solution.uniform_soln( self.mesh, 
                          node_elem_flag = 'ELEM', values = [2.0,'rand'] )

        ### Establish variables based on the mesh for testing later
        for soln in (self.node_solnA, self.node_solnB):
            soln.add_var( 'X', self.mesh.xyz_pt['X'] )
            soln.add_var( 'Y', self.mesh.xyz_pt['Y'] )
            soln.add_var( 'Z', self.mesh.xyz_pt['Z'] )
        for soln in (self.elem_solnA, self.elem_solnB):
            soln.add_var( 'X', self.mesh.xyz_el['X'] )
            soln.add_var( 'Y', self.mesh.xyz_el['Y'] )
            soln.add_var( 'Z', self.mesh.xyz_el['Z'] )

    def test_Max( self ):
        """Test the performance of soln_by_soln with operation="MAX"
        """

        ### Node-based solutions
        #   Test on var='q1'
        out_soln = steam.models.soln_by_soln( self.node_solnA, ['q1'], 
                                self.node_solnB, ['q1'], operation = 'MAX')
        self.assertTrue( all( out_soln.data['q1'] == 2 ) )
        #   Test max of X and Y in mesh
        self.assertTrue( all( self.node_solnA.data['X'] == 
                              self.mesh.xyz_pt['X'] ) )
        out_soln = steam.models.soln_by_soln( self.node_solnA, ['X'], 
                                self.node_solnB, ['Y'], operation = 'MAX')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_pt[['X','Y']].max(axis='columns') ) )
        #   Test max of X and Z in mesh
        self.assertTrue( all( self.node_solnA.data['X'] == 
                              self.mesh.xyz_pt['X'] ) )
        out_soln = steam.models.soln_by_soln( self.node_solnA, ['X'], 
                                self.node_solnB, ['Z'], operation = 'MAX')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_pt[['X','Z']].max(axis='columns') ) )

        ### Element-based solutions
        #   Test on var='q1'
        out_soln = steam.models.soln_by_soln( self.elem_solnA, ['q1'], 
                                self.elem_solnB, ['q1'], operation = 'MAX')
        self.assertTrue( all( out_soln.data['q1'] == 2 ) )
        #   Test max of X and Y in mesh
        self.assertTrue( all( self.elem_solnA.data['X'] == 
                              self.mesh.xyz_el['X'] ) )
        out_soln = steam.models.soln_by_soln( self.elem_solnA, ['X'], 
                                self.elem_solnB, ['Y'], operation = 'MAX')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_el[['X','Y']].max(axis='columns') ) )
        #   Test max of X and Z in mesh
        self.assertTrue( all( self.elem_solnA.data['X'] == 
                              self.mesh.xyz_el['X'] ) )
        out_soln = steam.models.soln_by_soln( self.elem_solnA, ['X'], 
                                self.elem_solnB, ['Z'], operation = 'MAX')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_el[['X','Z']].max(axis='columns') ) )

    def test_Min( self ):
        """Test the performance of soln_by_soln with operation="MIN"
        """

        ### Node-based solutions
        #   Test on var='q1'
        out_soln = steam.models.soln_by_soln( self.node_solnA, ['q1'], 
                                self.node_solnB, ['q1'], operation = 'MIN')
        self.assertTrue( all( out_soln.data['q1'] == 1 ) )
        #   Test min of X and Y in mesh
        self.assertTrue( all( self.node_solnA.data['X'] == 
                              self.mesh.xyz_pt['X'] ) )
        out_soln = steam.models.soln_by_soln( self.node_solnA, ['X'], 
                                self.node_solnB, ['Y'], operation = 'MIN')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_pt[['X','Y']].min(axis='columns') ) )
        #   Test min of X and Z in mesh
        self.assertTrue( all( self.node_solnA.data['X'] == 
                              self.mesh.xyz_pt['X'] ) )
        out_soln = steam.models.soln_by_soln( self.node_solnA, ['X'], 
                                self.node_solnB, ['Z'], operation = 'MIN')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_pt[['X','Z']].min(axis='columns') ) )

        ### Element-based solutions
        #   Test on var='q1'
        out_soln = steam.models.soln_by_soln( self.elem_solnA, ['q1'], 
                                self.elem_solnB, ['q1'], operation = 'MIN')
        self.assertTrue( all( out_soln.data['q1'] == 1 ) )
        #   Test min of X and Y in mesh
        self.assertTrue( all( self.elem_solnA.data['X'] == 
                              self.mesh.xyz_el['X'] ) )
        out_soln = steam.models.soln_by_soln( self.elem_solnA, ['X'], 
                                self.elem_solnB, ['Y'], operation = 'MIN')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_el[['X','Y']].min(axis='columns') ) )
        #   Test min of X and Z in mesh
        self.assertTrue( all( self.elem_solnA.data['X'] == 
                              self.mesh.xyz_el['X'] ) )
        out_soln = steam.models.soln_by_soln( self.elem_solnA, ['X'], 
                                self.elem_solnB, ['Z'], operation = 'MIN')
        self.assertTrue( all( out_soln.data.X == 
                         self.mesh.xyz_el[['X','Z']].min(axis='columns') ) )

    @ut.skip( 'Not implemented yet' )
    def test_Replacement( self ):
        """Test the functionality of the replacement model
        """

        ### Test replacement
        replace_factor = 4.5
        yfactors = [replace_factor] * len( self.y_gte0_nodes )
        steam.libmesh.equation_systems.model_at_nodes( self.es.solnp, 0, 
                            ModelType.REPLACEMENT, self.y_gte0_nodes, yfactors)
        
        xfactors = [replace_factor * 2] * len( self.x_gt0_nodes )
        steam.libmesh.equation_systems.model_at_nodes( self.es.solnp, 1, 
                            ModelType.REPLACEMENT, self.x_gt0_nodes, xfactors)
        
        ### Convert back to Python and probe for expected results
        out_soln = steam.solution.Solution.from_libmesh( self.es, 
                                                      name = self.sys_name )
        old_data_dict = self.soln.data.to_dict( orient = 'index' )
        new_data_dict = out_soln.data.to_dict( orient = 'index' )
        for node_num in new_data_dict.keys():
            ### Test replacement based on y
            if self.node_dict[node_num]['Y'] >= 0:
                self.assertTrue( new_data_dict[node_num]['q1'] == 
                                 replace_factor )
            else:
                self.assertTrue( new_data_dict[node_num]['q1'] == 
                                 old_data_dict[node_num]['q1'] )

            ### Test replacement based on x
            if self.node_dict[node_num]['X'] > 0:
                self.assertTrue( new_data_dict[node_num]['q2'] == 
                                 2 * replace_factor )
            else:
                self.assertTrue( new_data_dict[node_num]['q2'] == 
                                 old_data_dict[node_num]['q2'] )
        
    @ut.skip( 'Not implemented yet' )
    def test_Multiplicative( self ):
        """Test the multiplicative modeling functionality
        """

        ### Test multiplicative modeling
        mult_factor = 2.5
        yfactors = [mult_factor] * len( self.y_gte0_nodes )
        steam.libmesh.equation_systems.model_at_nodes( self.es.solnp, 0, 
                        ModelType.MULTIPLICATIVE, self.y_gte0_nodes, yfactors )
        
        xfactors = [mult_factor * 2] * len( self.x_gt0_nodes )
        steam.libmesh.equation_systems.model_at_nodes( self.es.solnp, 1, 
                        ModelType.MULTIPLICATIVE, self.x_gt0_nodes, xfactors )
        
        ### Write output
        out_file = 'multiplied.plt'
        self.es.write( 'PLT', out_file )
        self.assertTrue( os.path.isfile( out_file ) )
        self.remove_files.append( out_file )

        ### Convert back to Python and probe for expected results
        out_soln = steam.solution.Solution.from_libmesh( self.es, 
                                                      name = self.sys_name )
        old_data_dict = self.soln.data.to_dict( orient = 'index' )
        new_data_dict = out_soln.data.to_dict( orient = 'index' )
        for node_num in new_data_dict.keys():
            ### Test replacement based on y
            if self.node_dict[node_num]['Y'] >= 0:
                self.assertTrue( new_data_dict[node_num]['q1'] == 
                             mult_factor * old_data_dict[node_num]['q1'] )
            else:
                self.assertTrue( new_data_dict[node_num]['q1'] == 
                                 old_data_dict[node_num]['q1'] )

            ### Test replacement based on x
            if self.node_dict[node_num]['X'] > 0:
                self.assertTrue( new_data_dict[node_num]['q2'] == 
                             2 * mult_factor * old_data_dict[node_num]['q2'] )
            else:
                self.assertTrue( new_data_dict[node_num]['q2'] == 
                                 old_data_dict[node_num]['q2'] )
        
    @ut.skip( 'Not implemented yet' )
    def test_Additive( self ):
        """Test the additive modeling functionality
        """

        ### Test additive modeling
        add_factor = 5.0
        yfactors = [add_factor] * len( self.y_gte0_nodes )
        steam.libmesh.equation_systems.model_at_nodes( self.es.solnp, 0, 
                            ModelType.ADDITIVE, self.y_gte0_nodes, yfactors )
        
        xfactors = [add_factor * -1] * len( self.x_gt0_nodes )
        steam.libmesh.equation_systems.model_at_nodes( self.es.solnp, 1, 
                            ModelType.ADDITIVE, self.x_gt0_nodes, xfactors )
        
        ### Write output
        out_file = 'additive.plt'
        self.es.write( 'PLT', out_file )
        self.assertTrue( os.path.isfile( out_file ) )
        self.remove_files.append( out_file )

        ### Convert back to Python and probe for expected results
        out_soln = steam.solution.Solution.from_libmesh( self.es, 
                                                      name = self.sys_name )
        old_data_dict = self.soln.data.to_dict( orient = 'index' )
        new_data_dict = out_soln.data.to_dict( orient = 'index' )
        for node_num in new_data_dict.keys():
            ### Test replacement based on y
            if self.node_dict[node_num]['Y'] >= 0:
                self.assertTrue( new_data_dict[node_num]['q1'] == 
                                 old_data_dict[node_num]['q1'] + add_factor )
            else:
                self.assertTrue( new_data_dict[node_num]['q1'] == 
                                 old_data_dict[node_num]['q1'] )

            ### Test replacement based on x
            if self.node_dict[node_num]['X'] > 0:
                self.assertTrue( new_data_dict[node_num]['q2'] == 
                                 old_data_dict[node_num]['q2'] - add_factor )
            else:
                self.assertTrue( new_data_dict[node_num]['q2'] == 
                                 old_data_dict[node_num]['q2'] )
        
    def tearDown( self ):
        """Delete the mesh and solution object to ensure that independent
        tests don't affect each other.
        """ 

        self.mesh = None
        self.node_soln = None
        self.elem_soln = None
        self.model_soln = None

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """

        for f in cls.remove_files:
            if os.path.isfile( f ):
                os.remove( f )

#! Run the tests
if __name__ == '__main__':
    ut.main()

