#! /bin/env python

import unittest as ut
import steam

class TestProject( ut.TestCase ):


    def test_pickle( self ):
        """ Make a cube, add a solution, and just make sure that all of the
        I/O routines function.  Read-back any that we can and compare to the
        original.  Again, the goal is to avoid exceptions more than anything else.
        """
        mesh = steam.mesh.mk_cube(1.0,11)
        mesh.get_xyz_el()
        soln = steam.solution.uniform_soln( mesh, values = ['rand'])

        #! Mesh only
        steam.io.write_pickle(name="test.pkl",mesh=mesh,soln=None)
        (mesh2,soln2) = steam.io.read_pickle(file="test.pkl")

        self.assertTrue(mesh == mesh2)

        #! Both
        steam.io.write_pickle(name="test.pkl",mesh=mesh,soln=soln)
        (mesh2,soln2) = steam.io.read_pickle(file="test.pkl")

        self.assertTrue(mesh == mesh2)
        self.assertTrue(soln == soln2)
        pass


    def test_triq_ascii( self ):
        """ Make a cube, add a solution, and just make sure that all of the
        I/O routines function.  Read-back any that we can and compare to the
        original.  Again, the goal is to avoid exceptions more than anything else.

        This is triq, so node-based.
        """
        return
        mesh = steam.mesh.mk_cube(1.0,11)
        soln = steam.solution.uniform_soln( mesh, "NODE",values = ['rand'])

        #! Mesh only
        steam.io.write_triq_ascii(name="test.triq",mesh=mesh,soln=None)
        (mesh2,soln2) = steam.io.read_triq_ascii(file="test.triq",toElement=False)

        self.assertTrue(mesh == mesh2)

        #! Both
        steam.io.write_triq_ascii(name="test.triq",mesh=mesh,soln=soln)
        (mesh2,soln2) = steam.io.read_triq_ascii(file="test.triq",toElement=False)

        self.assertTrue(mesh == mesh2)
        self.assertTrue(soln == soln2)
        pass

    def test_triq_uform( self ):
        """ Make a cube, add a solution, and just make sure that all of the
        I/O routines function.  Read-back any that we can and compare to the
        original.  Again, the goal is to avoid exceptions more than anything else.

        This is triq, so node-based.
        """
        mesh = steam.mesh.mk_cube(1.0,11)
        soln = steam.solution.uniform_soln( mesh, "NODE",values = ['rand'])

        #! Mesh only
        steam.io.write_triq_uform(name="test.triq",mesh=mesh,soln=None)
        (mesh2,soln2) = steam.io.read_triq_uform(file="test.triq",toElement=False)

        self.assertTrue(mesh == mesh2)

        #! Both
        steam.io.write_triq_uform(name="test.triq",mesh=mesh,soln=soln)
        (mesh2,soln2) = steam.io.read_triq_uform(file="test.triq",toElement=False)

        self.assertTrue(mesh == mesh2)


        #! Single precision should NOT be identical

        #! Mesh only
        steam.io.write_triq_uform(name="test.triq",mesh=mesh,soln=None,dp=False)
        (mesh2,soln2) = steam.io.read_triq_uform(file="test.triq",toElement=False,dp=False)

        self.assertFalse(mesh == mesh2)

        #! Both
        steam.io.write_triq_uform(name="test.triq",mesh=mesh,soln=soln,dp=False)
        (mesh2,soln2) = steam.io.read_triq_uform(file="test.triq",toElement=False,dp=False)

        self.assertFalse(mesh == mesh2)
        self.assertFalse(soln == soln2)


        pass

    def test_vtu( self ):
        """ Make a cube, add a solution, and just make sure that all of the
        I/O routines function.  The goal is to avoid exceptions more than anything else.
        """

        #Node Based
        mesh = steam.mesh.mk_cube(1.0,11)
        soln = steam.solution.uniform_soln( mesh, "NODE",values = ['rand'])

        steam.io.write_vtu(name="test.vtu",mesh=mesh,soln=None)
        steam.io.write_vtu(name="test.vtu",mesh=mesh,soln=soln)
        steam.io.write_vtu(name="test.vtu",soln=soln)

        #Element Based
        mesh = steam.mesh.mk_cube(1.0,11)
        soln = steam.solution.uniform_soln( mesh,values = ['rand'])

        steam.io.write_vtu(name="test.vtu",mesh=mesh,soln=None)
        steam.io.write_vtu(name="test.vtu",mesh=mesh,soln=soln)
        steam.io.write_vtu(name="test.vtu",soln=soln)

        pass

    def test_cdat( self ):
        """ Make a cube, add a solution, and just make sure that all of the
        I/O routines function.  The goal is to avoid exceptions more than anything else.
        """
        mesh = steam.mesh.mk_cube(1.0,11)
        soln = steam.solution.uniform_soln( mesh, "NODE",values = ['rand'])

        #! Mesh only
        steam.io.write_vtu(name="test.cdat",mesh=mesh,soln=None)
        steam.io.write_vtu(name="test.cdat",mesh=mesh,soln=soln)
        steam.io.write_vtu(name="test.cdat",soln=soln)
        pass

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

        #files_to_remove = ( 'test.dat', 'test.e', 'original_soln.triq', 
        #                    'new_soln.triq', 'cube.tri' )
        #for f in files_to_remove:
        #    if os.path.isfile( f ):
        #        os.remove( f )


#! Run the tests
if __name__ == '__main__':
    ut.main()

