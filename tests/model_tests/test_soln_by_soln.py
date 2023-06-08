#! /bin/env python

""" Replicate the LibMESH tests for the Python-side modeling methods.
"""

import unittest as ut
import steam
import os
from code import interact

@ut.skip( 'Not implemented yet' )
class TestModels( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        #self.mesh.write_tri_ascii( 'cube' )
        cls.remove_files = []

        ### Initialize LibMesh
        steam.libmesh.Init()

    def setUp( self ):
        """Create a 5x5x5 cube mesh with two uniform solution variables.
        """

        ### Create a cube and populate it with uniform solutions
        self.mesh = steam.mesh.mk_cube( 1, 5, file = 0 )
        self.node_soln = steam.solution.uniform_soln( self.mesh, 
                           node_elem_flag = 'NODE', values = ['rand','rand'] )
        self.elem_soln = steam.solution.uniform_soln( self.mesh, 
                           node_elem_flag = 'ELEM', values = ['rand','rand'] )

        ### Define a dictionary of all node_ids and their associated points
        self.node_dict = self.mesh.xyz_pt.to_dict( orient = 'index' )
        self.node_ids = list( map( int, self.node_dict.keys() ) ) 

        ### Define lists of points with various qualities
        self.x_gt0_nodes = list( map( int, [node for node in 
                self.node_dict.keys() if self.node_dict[node]['X'] > 0] ) )
        self.y_gte0_nodes = list( map( int, [node for node in 
                self.node_dict.keys() if self.node_dict[node]['Y'] >= 0] ) )

        ### Construct the modeling solution
        self.elem_model_soln = steam.solution.uniform_soln( self.mesh, 
                       node_elem_flag = 'ELEM', values = [1.0, 1.0, 0.0] )
        self.elem_model_soln.data.rename( columns = {'q1': 'repl', 
                                          'q2': 'mult', 'q3': 'add' } )
        steam.models.scalar_to_component( self.elem_model_soln, 1, 'repl',
                                          self.replace_factor, method='repl' )

        interact( local = dict( globals(), **locals() ) )

    def test_Setup( self ):
        """Confirm that the setUp method was run properly.
        """
    
        out_file = 'initial_soln.plt'

        self.assertEqual( len( self.x_gt0_nodes ), 41 )
        self.assertEqual( len( self.y_gte0_nodes ), 57 )

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

