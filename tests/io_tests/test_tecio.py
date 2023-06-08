#! /bin/env python

import unittest as ut
import steam
has_pydata = False
try:
    import pydata
    has_pydata = True
except ImportError:
    pass

class TestProject( ut.TestCase ):


    @ut.skipUnless( has_pydata, 'pydata is not present in PYTHONPATH' )
    def test_tecplot_write_grid( self ):
        """ Make a cube, write it to tecplot format, then read it back in and verify it is the same as the original.
        """
        #import pydata
        import numpy as np


        tolerance = 1.0e-10
        all_nodes_equal = True
        all_elems_equal = True

        mesh = steam.mesh.mk_cube(1.0,11)

        #! Mesh only
        steam.io.write_tecplot_ascii(name="test.dat",mesh=mesh)
        (mesh2,soln2) = steam.io.read_tecplot_ascii(file="test.dat")

        map = {}
        # Check xyz_pt
        for i in range (len (mesh.xyz_pt.values.tolist ())):
            node1 = mesh.xyz_pt.values.tolist ()[i]
            j = -1
            node_equal = False
            while not node_equal and j < len (mesh2.xyz_pt.values.tolist ()) - 1:
                j += 1
                node2 = mesh2.xyz_pt.values.tolist ()[j]
                node_equal = True
                for n in range (len (node1)):
                    if abs (node1 [n] - node2 [n]) > tolerance: node_equal = False
                if node_equal: map [j] = i
            if not node_equal:
                all_nodes_equal = False
                print ("Node {} - {} did not find a matching node within tolerance ({})".format (i, node1, tolerance))

#        for key, value in map.items():
#            if key != value: print (key, value)

#        Did not work since nodes are not ordered the same way
#        self.assertTrue(np.allclose(mesh.xyz_pt.values, mesh2.xyz_pt.values, rtol=1.0e-10, atol=1.0e-8))

        self.assertTrue(all_nodes_equal)

        # Check conn
        for i in range (len (mesh.conn.values.tolist ())):
            elem = mesh.conn.values.tolist ()[i]
            j = -1
            elem_equal = False
            while not elem_equal and j < len (mesh2.conn.values.tolist ()) - 1:
                j += 1
                elem2 = mesh2.conn.values.tolist ()[j]
                elem_equal = True
                for n in range (len (elem2)):
                    if elem [n] != map [elem2 [n]]: elem_equal = False
            if not elem_equal:
                all_elems_equal = False
                print ("Element {} - {} did not find a match".format (i, elem))

        self.assertTrue(all_elems_equal)

        pass

    @ut.skipUnless( has_pydata, 'pydata is not present in PYTHONPATH' )
    def test_tecplot_read_write_point_cloud( self ):
        """ Make a cube, set it as point_cloud write it to tecplot format, 
            then read it back in and verify it is the same as the original.
        """
        #import pydata
        import numpy as np

        tolerance = 1.0e-10
        all_nodes_equal = True
        all_elems_equal = True

        mesh = steam.mesh.mk_cube(1.0,11)
        mesh.to_point_cloud( keep_nodes = True )

        #! Mesh only
        steam.io.write_tecplot_ascii(name="point_cloud.dat",mesh=mesh)
        (mesh2,soln2) = steam.io.read_tecplot_ascii(file="point_cloud.dat",pointcloud=True)

        self.assertTrue(np.allclose(mesh.xyz_pt.values, mesh2.xyz_pt.values, rtol=1.0e-10, atol=1.0e-8))


    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """
        pass
        import os

        #os.remove('file')
        #import shutil
        #shutil.rmtree('dir')

        files_to_remove = ( 'test.dat', 'point_cloud.dat' )
        for f in files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )

#! Run the tests
if __name__ == '__main__':
    ut.main()

