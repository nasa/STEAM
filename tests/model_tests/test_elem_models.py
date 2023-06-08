#! /bin/env python

import unittest as ut
import steam
import os
if steam.has_libmesh:
    from steam.libmesh.equation_systems import ModelType
else:
    ModelType = None

@ut.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
class TestProject( ut.TestCase ):

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
        self.soln = steam.solution.uniform_soln( self.mesh, 
                           node_elem_flag = 'ELEM', values = ['rand','rand'] )

        ### Put it in libMesh
        self.lmesh = steam.libmesh.Mesh.from_python( self.mesh )
        self.es = steam.libmesh.EquationSystems( self.lmesh )
        self.sys_name = self.es.from_python( self.soln )

        #### Define a dictionary of all element centroids
        self.mesh.get_xyz_el()
        self.centroid_dict = self.mesh.xyz_el.to_dict( orient = 'index' )

        ### Define lists of elements with various qualities
        self.x_gt0_elems = list( map( int, [elem for elem in 
                self.centroid_dict.keys() 
                if self.centroid_dict[elem]['X'] > 0] ) )
        self.y_gte0_elems = list( map( int, [elem for elem in 
                self.centroid_dict.keys() 
                if self.centroid_dict[elem]['Y'] >= 0] ) )

    def test_Setup( self ):
        """Confirm that the setUp method was run properly.
        """
    
        out_file = 'initial_soln.plt'
        self.es.write( 'PLT', out_file )
        self.assertTrue( os.path.isfile( out_file ) )
        self.remove_files.append( out_file )

        self.assertEqual( len( self.x_gt0_elems ), 96 )
        self.assertEqual( len( self.y_gte0_elems ), 96 )

    def test_Replacement( self ):
        """Test the functionality of the replacement model
        """

        ### Test replacement
        replace_factor = 4.5
        yfactors = [replace_factor] * len( self.y_gte0_elems )
        steam.libmesh.equation_systems.model_at_elems( self.es.solnp, 0, 
                        ModelType.REPLACEMENT, self.y_gte0_elems, yfactors )
        xfactors = [replace_factor * 2] * len( self.x_gt0_elems )
        steam.libmesh.equation_systems.model_at_elems( self.es.solnp, 1, 
                        ModelType.REPLACEMENT, self.x_gt0_elems, xfactors )
        
        ### Write output
        out_file = 'replaced.plt'
        self.es.write( 'PLT', out_file )
        self.assertTrue( os.path.isfile( out_file ) )
        self.remove_files.append( out_file )

        ### Convert back to Python and probe for expected results
        out_soln = steam.solution.Solution.from_libmesh( self.es, 
                                                      name = self.sys_name )
        old_data_dict = self.soln.data.to_dict( orient = 'index' )
        new_data_dict = out_soln.data.to_dict( orient = 'index' )
        for elem_num in new_data_dict.keys():
            ### Test replacement based on y
            if self.centroid_dict[elem_num]['Y'] >= 0:
                self.assertTrue( new_data_dict[elem_num]['q1'] == 
                                 replace_factor )
            else:
                self.assertTrue( new_data_dict[elem_num]['q1'] == 
                                 old_data_dict[elem_num]['q1'] )

            ### Test replacement based on x
            if self.centroid_dict[elem_num]['X'] > 0:
                self.assertTrue( new_data_dict[elem_num]['q2'] == 
                                 2 * replace_factor )
            else:
                self.assertTrue( new_data_dict[elem_num]['q2'] == 
                                 old_data_dict[elem_num]['q2'] )
        
    def test_Multiplicative( self ):
        """Test the multiplicative modeling functionality
        """

        ### Test replacement
        mult_factor = 2.0
        yfactors = [mult_factor] * len( self.y_gte0_elems )
        steam.libmesh.equation_systems.model_at_elems( self.es.solnp, 0, 
                    ModelType.MULTIPLICATIVE, self.y_gte0_elems, yfactors )
        
        xfactors = [mult_factor * 2] * len( self.x_gt0_elems )
        steam.libmesh.equation_systems.model_at_elems( self.es.solnp, 1, 
                    ModelType.MULTIPLICATIVE, self.x_gt0_elems, xfactors )
        
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
        for elem_num in new_data_dict.keys():
            ### Test replacement based on y
            if self.centroid_dict[elem_num]['Y'] >= 0:
                self.assertTrue( new_data_dict[elem_num]['q1'] == 
                             mult_factor * old_data_dict[elem_num]['q1'] )
            else:
                self.assertTrue( new_data_dict[elem_num]['q1'] == 
                             old_data_dict[elem_num]['q1'] )

            ### Test replacement based on x
            if self.centroid_dict[elem_num]['X'] > 0:
                self.assertTrue( new_data_dict[elem_num]['q2'] == 
                             2 * mult_factor * old_data_dict[elem_num]['q2'] )
            else:
                self.assertTrue( new_data_dict[elem_num]['q2'] == 
                             old_data_dict[elem_num]['q2'] )
        
    def test_Additive( self ):
        """Test the additive modeling functionality
        """

        ### Test replacement
        add_factor = 5.0
        yfactors = [add_factor] * len( self.y_gte0_elems )
        steam.libmesh.equation_systems.model_at_elems( self.es.solnp, 0, 
                            ModelType.ADDITIVE, self.y_gte0_elems, yfactors )
        
        xfactors = [add_factor * -1] * len( self.x_gt0_elems )
        steam.libmesh.equation_systems.model_at_elems( self.es.solnp, 1, 
                            ModelType.ADDITIVE, self.x_gt0_elems, xfactors )
        
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
        for elem_num in new_data_dict.keys():
            ### Test addition based on y
            if self.centroid_dict[elem_num]['Y'] >= 0:
                self.assertTrue( new_data_dict[elem_num]['q1'] == 
                                 old_data_dict[elem_num]['q1'] + add_factor )
            else:
                self.assertTrue( new_data_dict[elem_num]['q1'] == 
                                 old_data_dict[elem_num]['q1'] )

            ### Test addition based on x
            if self.centroid_dict[elem_num]['X'] > 0:
                self.assertTrue( new_data_dict[elem_num]['q2'] == 
                                 old_data_dict[elem_num]['q2'] - add_factor )
            else:
                self.assertTrue( new_data_dict[elem_num]['q2'] == 
                                 old_data_dict[elem_num]['q2'] )
        
    def tearDown( self ):
        """Delete the mesh and solution object to ensure that independent
        tests don't affect each other.
        """ 

        self.mesh = None
        self.soln = None
        self.lmesh = None
        self.es = None

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

