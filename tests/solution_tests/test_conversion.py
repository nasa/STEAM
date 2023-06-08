#! /bin/env python

import unittest as ut
import steam
import pandas as pd

class TestProject( ut.TestCase ):


    @classmethod
    def setUpClass( cls ):
        """Create temporary files. This method gets called once, at the very
        beginning of testing.
        """

        #! Make a simple tri mesh and solution with fake data
        cls.mesh = steam.mesh.Mesh()
        nodes = {'X' : [0, 1, 2, 3], 'Y' : [0, 1, 1, 2], 'Z' : [0, 0, 0, 0]}
        elems = {'0' : [0, 1],
                 '1' : [1, 2],
                 '2' : [2, 3]}
        comps = {'0' : [1],
                 '1' : [1]}

        cls.mesh.xyz_pt =  pd.DataFrame.from_dict (nodes)
        cls.mesh.conn   = pd.DataFrame.from_dict (elems)
        cls.mesh.comp   = pd.DataFrame.from_dict (comps, 'index')

        cls.nsoln =  steam.solution.Solution(mesh=cls.mesh)
        cls.nsoln.data = pd.DataFrame.from_dict (nodes)

        cls.esoln = steam.solution.Solution(mesh=cls.mesh)
        cls.esoln.data = pd.DataFrame.from_dict ({'X' : [1.0    , 2.0], 
                                                  'Y' : [2.0/3.0, 4.0/3.0], 
                                                  'Z' : [0.0    , 0.0]})

        #! Ok, I did the math by hand here for element to node
        cls.e2n = pd.DataFrame.from_dict ({'X' : [1.0,     1.5,     1.5 ,    2.0    ], 
                                           'Y' : [2.0/3.0, 3.0/3.0, 3.0/3.0, 4.0/3.0], 
                                           'Z' : [0.0    , 0.0,     0.0,     0.0    ]
                                          })
        cls.nsoln.point   = True
        cls.esoln.element = True

        pass

    def test_write_form( self ):
        """ Test formated read/write of triq file.
        """

        #! Mapping of variable names
        mvars = ['X','Y','Z']
        fvars = ['q1','q2','q3']

        #! Write the node-based data out
        steam.io.write_triq_ascii(name="node.triq",soln=self.nsoln)
        (tmesh, tmp) = steam.io.read_triq_ascii("node.triq",toElement=False)

        for (i,row) in tmp.data.iterrows():
            for v in range(3):
                mv = mvars[v]
                fv = fvars[v]
                self.assertAlmostEqual(row[fv],self.nsoln.data[mv][i])


        #! Write the element-based data out (converts it to nodes
        steam.io.write_triq_ascii(name="elem.triq",soln=self.esoln)
        (tmesh, tmp) = steam.io.read_triq_ascii("elem.triq",toElement=False)

        for (i,row) in self.e2n.iterrows():
            for v in range(3):
                mv = mvars[v]
                fv = fvars[v]
                self.assertAlmostEqual(row[mv],tmp.data[fv][i])

    def test_write_uform( self ):
        """ Test unformated read/write of triq file.
        """

        #! Mapping of variable names
        mvars = ['X','Y','Z']
        fvars = ['q1','q2','q3']

        #! Write the node-based data out
        steam.io.write_triq_uform(name="node2.triq",soln=self.nsoln)
        (tmesh, tmp) = steam.io.read_triq_uform("node2.triq",toElement=False)

        for (i,row) in tmp.data.iterrows():
            for v in range(3):
                mv = mvars[v]
                fv = fvars[v]
                self.assertAlmostEqual(row[fv],self.nsoln.data[mv][i])


        #! Write the element-based data out (converts it to nodes
        steam.io.write_triq_uform(name="elem2.triq",soln=self.esoln)
        (tmesh, tmp) = steam.io.read_triq_uform("elem2.triq",toElement=False)

        for (i,row) in self.e2n.iterrows():
            for v in range(3):
                mv = mvars[v]
                fv = fvars[v]
                self.assertAlmostEqual(row[mv],tmp.data[fv][i])


    def test_to_elem( self ):
        """ Convert from node to elements.

        """

        tmp = self.nsoln.node_to_element(return_new=True)

        for (i,row) in tmp.data.iterrows():
            for var in ['X','Y','Z']:
        #       print("Comparing ",i," ",var)
                self.assertAlmostEqual(row[var],self.esoln.data[var][i])


    def test_to_node( self ):
        """ Convert from elements to nodes.

        """

        tmp = self.esoln.element_to_node(return_new=True)

        for (i,row) in self.e2n.iterrows():
            for var in ['X','Y','Z']:
        #       print("Comparing ",i," ",var)
                self.assertAlmostEqual(row[var],tmp.data[var][i])

    @ut.skipUnless(steam.has_libmesh, "test requires a compiled LibMesh")
    def test_python_to_libmesh(self):
        """ Test the python-to-libmesh and libmesh-to-python conversion.  

        NOTE: This test is reproduced by ../libmesh_tests/test_solution.py
        """

        ### Initialize LibMesh
        steam.libmesh.Init()

        ### Make the cube and convert it to libmesh
        cube = steam.mesh.mk_cube(5.0,31)
        lcube = steam.libmesh.Mesh.from_python(cube)

        ### Make the solution and convert it to libmesh
        soln = steam.solution.uniform_soln( cube,
                                   node_elem_flag = 'ELEM', values = [1,3] )
        es = steam.libmesh.EquationSystems( lcube )
        sys_name = es.from_python( soln )


        ### Convert it back
        #cube2 = steam.mesh.Mesh.from_libmesh(lcube)
        soln2 = steam.solution.Solution.from_libmesh( es, name = sys_name )

        #### Compare everything
        ## Points are the same:
        ## - first  all() collapses rows
        ## - second all() collapses columns
        #all_points = ((cube.points == cube2.points).all()).all()
        #self.assertTrue(all_points,"Points do not match")
        ## Connectivity is  the same:
        #all_conn   = ((cube.ele    == cube2.ele   ).all()).all()
        #self.assertTrue(all_conn  ,"Connectivity does not match")
        ## Components   are the same:
        #all_comp   = ((cube.comp   == cube2.comp  ).all()).all()
        #self.assertTrue(all_comp  ,"Components do not match")
        
        all_data = ((soln.data == soln2.data).all()).all()
        self.assertTrue(all_data,"Data does not match")
        self.assertTrue(soln.data.equals(soln2.data),"Data does not match")

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """
        pass
        #import os
        #os.remove('file')
        #import shutil
        #shutil.rmtree('dir')


#! Run the tests
if __name__ == '__main__':
    ut.main()

