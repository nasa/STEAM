#! /bin/env python

import unittest as ut
import steam
import numpy as np

class TestProject( ut.TestCase ):

    def test_read_write( self):
        """ Make a transformation, use it, and then read/write it and use it again."""

        mesh = steam.mesh.mk_cube( 1, 15)
        mesh.get_xyz_el()
        soln = steam.solution.uniform_soln( mesh, node_elem_flag = 'NODE',
                                            values = [1.0,10.0,20.0])
        soln.node_to_element()
        
        steam.models.scalar_to_component(soln,1,'q1',1.0)
        steam.models.scalar_to_component(soln,2,'q1',2.0)
        steam.models.scalar_to_component(soln,3,'q1',3.0)
        steam.models.scalar_to_component(soln,4,'q1',4.0,'repl')
        steam.models.scalar_to_component(soln,5,'q1',5.0,'repl')
        steam.models.scalar_to_component(soln,6,'q1',6.0,'repl')
        
        target = steam.mesh.mk_sphere( .8, 25 )
        target.get_xyz_el()
        trans =  steam.interpolate.Transform(mesh,target)
        dist  =  trans.inverse_distance(k=5,p=2)
        
        trans_quad =  steam.interpolate.Transform(mesh,target)
        dist       =  trans_quad.inverse_distance_quad(k=5,p=2)
        
        mem_soln      = trans.apply(soln)
        mem_soln_quad = trans_quad.apply(soln)
        
        cont = steam.container.Container("test.h5","w")
        cont.add("/trans",trans)
        cont.add("/trans_quad",trans_quad)
        cont.write()
        
        cont = steam.container.Container("test.h5")
        cont.read()
        
        disk_soln      = cont.obj['trans'].apply(soln)
        disk_soln_quad = cont.obj['trans_quad'].apply(soln)

        #! The ones in memory should agree to the ones from the disk transform
        check =np.all(mem_soln.data.values == disk_soln.data.values)
        self.assertTrue( check )
        check =np.all(mem_soln_quad.data.values == disk_soln_quad.data.values)
        self.assertTrue( check )

        #! Quad should not equal non-quad
        check =np.all(mem_soln_quad.data.values == disk_soln.data.values)
        self.assertFalse( check )
        check =np.all(mem_soln.data.values == disk_soln_quad.data.values)
        self.assertFalse( check )

        return

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """
        import os

        files_to_remove = ( 'test.h5')
        for f in files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )
#! Run the tests
if __name__ == '__main__':
    ut.main()

