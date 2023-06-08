#! /bin/env python

import unittest as ut
import steam
import numpy as np
import pandas as pd

def get_simple_cube(dim,l=1.0):
    """Make a cube with Pressure on the min sides."""
    mesh = steam.mesh.mk_cube( l, dim)
    mesh.get_xyz_el()
    soln = steam.solution.uniform_soln( mesh, values = [1.0])
    soln.rename_vars({'q1':'cp'})
    steam.models.scalar_to_component(soln,1,'cp',1.0)
    steam.models.scalar_to_component(soln,2,'cp',0.0)
    steam.models.scalar_to_component(soln,3,'cp',2.0)
    steam.models.scalar_to_component(soln,4,'cp',0.0)
    steam.models.scalar_to_component(soln,5,'cp',3.0)
    steam.models.scalar_to_component(soln,6,'cp',0.0)
    return (mesh,soln)

class TestFomoco( ut.TestCase ):

    def test_analytical( self ):
        """ This test will to the math by hand and compare.
        """

        # Just integrate the faces
        (mesh,soln) = get_simple_cube(11)
        (force,moment) = steam.aero_util.fomoco(soln)
        np.testing.assert_allclose(force ,[1.0,2.0,3.0])
        np.testing.assert_allclose(moment,0,atol=1e-16)
    
        # So, the absolute values for forces worked out (above)
        # Now more the MRC around
        weights = steam.aero_util.fomoco_weights(mesh)

        # Just integrate the faces
        (mesh,soln) = get_simple_cube(11)

        for xloc in [-0.5,0.0,-1,2,20]:
            for yloc in [-0.5,0.0,-1,2,20]:
                for zloc in [-0.5,0.0,-1,2,20]:
                    (force,moment) = steam.aero_util.fomoco(soln,mrc=[xloc,yloc,zloc],weights=weights)
                    np.testing.assert_allclose(force ,[1.0,2.0,3.0])
                    np.testing.assert_allclose(moment,
                                   [zloc*force[1]-yloc*force[2],
                                    xloc*force[2]-zloc*force[0],
                                    yloc*force[0]-xloc*force[1]
                                    ],
                                   atol=1e-16
                                   )

        for xloc in np.random.rand(3)-.5:
            for yloc in np.random.rand(3)-.5:
                for zloc in np.random.rand(3)-.5:
                    (force,moment) = steam.aero_util.fomoco(soln,mrc=[xloc,yloc,zloc],weights=weights)
                    np.testing.assert_allclose(force ,[1.0,2.0,3.0])
                    np.testing.assert_allclose(moment,
                                   [zloc*force[1]-yloc*force[2],
                                    xloc*force[2]-zloc*force[0],
                                    yloc*force[0]-xloc*force[1]
                                    ],
                                   atol=1e-16
                                   )
        return

#! Run the tests
if __name__ == '__main__':
    ut.main()
