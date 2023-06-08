#! /bin/env python

import unittest as ut
import numpy as np
import steam

class TestLineLoads( ut.TestCase ):

    def test_check_calc( self ):
        """ Slice and check the returned ids."""
        
        mesh = steam.mesh.mk_cube( 1, 5)
        mesh.get_xyz_el()
        soln = steam.solution.uniform_soln(mesh,'ELEMS',values=['rand'])
        soln.rename_vars({'q1':'cp'})

        As = [1.0]
        Ls = [1.0]

        for dir  in ['X','Y','Z']:
            for (A,L)  in zip(As,Ls):
                for mrc in [[0,0,0],[.2,.5,1.0]]:
                    for res in [.1,1.0,0.5]:
                        ll = steam.aero_util.LineLoads(mesh,
                                                       A=A,
                                                       L=L,
                                                       mrc=mrc,
                                                       res=res,
                                                       dir=dir)
                        self.assertTrue(ll.check(soln))

    def test_write_ll (self):
        """ Write a line loads file, then read it and make sure it was lossless."""
        
        mesh = steam.mesh.mk_cube( 1, 5)
        mesh.get_xyz_el()
        soln = steam.solution.uniform_soln(mesh,'ELEMS',values=['rand'])
        soln.rename_vars({'q1':'cp'})
        ll = steam.aero_util.LineLoads(mesh,res=.05)
        fm_mem = ll.calc(soln)
        ll.write()
        fm_disk = np.loadtxt('output.ll')

        # OK, now compare to what's in memory first:
        np.testing.assert_allclose(fm_mem,ll.fomoco)

        # Now compare to what came from disk
        np.testing.assert_allclose(fm_mem     ,fm_disk[:,1:])
        np.testing.assert_allclose(ll.fomoco  ,fm_disk[:,1:])
        np.testing.assert_allclose(ll.centroids,fm_disk[:,0])


pass
        
        
#! Run the tests
if __name__ == '__main__':
    ut.main()

