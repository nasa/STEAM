#! /bin/env python

import unittest as ut
import numpy as np
import steam

class TestSlice( ut.TestCase ):

    def test_tri_count( self ):
        """ Slice between nodes to make new tris. Don't touch components"""
        ### Make the cubes
        orig = steam.mesh.mk_cube( 1, 5)

        mesh = steam.mesh.mk_cube( 1, 5)
        # One cut, on a seam - no new tris
        planes = np.zeros((1,3,3))
        planes[0,1,1] = 1
        planes[0,2,2] = 1
        #idlist = mesh.slice([[[0,0,0],[0,1,0],[0,0,1]]])
        idlist = mesh.slice(planes)

        # This should cut nicely around the cube, so we'll cut
        # 8 tris on each face.  Each cut tri becomes 3 tris.
        # 8*3*4 = 96 new triangles, but 32 tris go away.
        self.assertEqual(
                        orig.conn.values.shape[0],
                        mesh.conn.values.shape[0]
                        )

        # Ok, make some new tris
        #idlist = mesh.slice([[[.1,0,0],[.1,1,0],[.1,0,1]]])
        planes[:,:,0] = .1
        idlist = mesh.slice(planes)

        # This should cut nicely around the cube, so we'll cut
        # 8 tris on each face.  Each cut tri becomes 3 tris.
        # 8*3*4 = 96 new triangles, but 32 tris go away.
        self.assertEqual(
                        orig.conn.values.shape[0]+ 64,
                        mesh.conn.values.shape[0]
                        )

        # Do it again with two cuts
        mesh = steam.mesh.mk_cube( 1, 5)
        idlist = mesh.slice([[[.1,0,0],[.1,1,0],[.1,0,1]],
                             [[0,.2,0],[1,.2,0],[0,.2,1]],
            ])

        # This should cut nicely around the cube, so we'll cut
        # 8 tris on each face.  Each cut tri becomes 3 tris.
        # 8*3*4 = 96 new triangles, but 32 tris go away.
        # The second cut, does the same thing, but now there
        # are 11*2*3+8*3*2 new tris (minus 11*2+8*2 original)
        self.assertEqual(
                        orig.conn.values.shape[0]+ 64 + 44 + 32,
                        mesh.conn.values.shape[0]
                        )

        # Do it again with two cuts, but hit a vertex
        mesh = steam.mesh.mk_cube( 1, 5)
        idlist = mesh.slice([[[.1,0,0],[.1,1,0],[.1,0,1]],
                             [[0,.1,0],[1,.1,0],[0,.1,1]],
            ])

        self.assertEqual(
                        orig.conn.values.shape[0]+ 132,
                        mesh.conn.values.shape[0]
                        )

        #steam.io.write_vtu(name="test.vtu",mesh=mesh)

    def test_delete( self ):
        """ Slice and delete portions of the mesh"""

        ### Make the original for comparison
        orig = steam.mesh.mk_cube( 1, 5)

        # Make a cut, and delete everything on one side
        # of it.  One plane.  This does not create any new triangles
        # since we are slicing on an edge:
        xplane = np.zeros((1,3,3))
        xplane[0,1,1] = 1
        xplane[0,2,2] = 1
        yplane = np.zeros((1,3,3))
        yplane[0,0,0] = 1
        yplane[0,2,2] = 1
        zplane = np.zeros((1,3,3))
        zplane[0,0,0] = 1
        zplane[0,1,1] = 1

        # Cut off the tris on either side of 0
        for (i,plane) in enumerate([xplane,yplane,zplane]):
            mesh = steam.mesh.mk_cube( 1, 5)
            idlist = mesh.slice(plane,delete=[(-1,1)])

            self.assertFalse(np.any(mesh.xyz_pt.values[:,i] < 0))
            self.assertEqual(
                            orig.conn.values.shape[0]/2,
                            mesh.conn.values.shape[0]
                            )

            mesh = steam.mesh.mk_cube( 1, 5)
            idlist = mesh.slice(plane,delete=[(1,-1)])

            self.assertFalse(np.any(mesh.xyz_pt.values[:,i] > 0))
            self.assertEqual(
                            orig.conn.values.shape[0]/2,
                            mesh.conn.values.shape[0]
                            )

        # OK, now actually cut thing not at a seam
        for (i,plane) in enumerate([xplane,yplane,zplane]):
            offset = .375
            plane[0,:,i] = offset
            mesh = steam.mesh.mk_cube( 1, 5)
            idlist = mesh.slice(plane,delete=[(-1,1)])

            steam.io.write_vtu(name="test.vtu",mesh=mesh)
            self.assertFalse(np.any(mesh.xyz_pt.values[:,i] < offset))

            mesh = steam.mesh.mk_cube( 1, 5)
            idlist = mesh.slice(plane,delete=[(1,-1)])

            self.assertFalse(np.any(mesh.xyz_pt.values[:,i] > offset))
                
    def test_del_with_el( self ):
        """ Slice and delete with xyz_el populated """

        ### Make the original for comparison
        orig = steam.mesh.mk_sphere( 1, 11)

        # Make a cut, and delete everything on one side
        # of it.  One plane.  This does not create any new triangles
        # since we are slicing on an edge:
        xplane = np.zeros((1,3,3))
        xplane[0,1,1] = 1
        xplane[0,2,2] = 1
        yplane = np.zeros((1,3,3))
        yplane[0,0,0] = 1
        yplane[0,2,2] = 1
        zplane = np.zeros((1,3,3))
        zplane[0,0,0] = 1
        zplane[0,1,1] = 1

        # Cut off the tris on either side of 0
        for (i,plane) in enumerate([xplane,yplane,zplane]):
            mesh = steam.mesh.Mesh( copy=orig )
            mesh.get_xyz_el()
            idlist = mesh.slice(plane,delete=[(-1,1)])

#            self.assertFalse(np.any(mesh.xyz_pt.values[:,i] < 0))
#            self.assertEqual(
#                            orig.conn.values.shape[0]/2,
#                            mesh.conn.values.shape[0]
#                            )
#
#            mesh = steam.mesh.mk_cube( 1, 5)
#            idlist = mesh.slice(plane,delete=[(1,-1)])
#
#            self.assertFalse(np.any(mesh.xyz_pt.values[:,i] > 0))
#            self.assertEqual(
#                            orig.conn.values.shape[0]/2,
#                            mesh.conn.values.shape[0]
#                            )


    def test_component( self ):
        """ Slice and change the components on either side"""
        plane = np.zeros((1,3,3))
        plane[0,1,1] = 1
        plane[0,2,2] = 1
        plane[0,:,0] = .1
        
        mesh = steam.mesh.mk_cube( 1, 5)
        mesh.get_xyz_el()
        idlist = mesh.slice(plane,component=[(-1,101)])
        np.testing.assert_equal(
                                mesh.comp.values[mesh.xyz_el.values[:,0] <.1],
                                101
                                )
        
        mesh = steam.mesh.mk_cube( 1, 5)
        mesh.get_xyz_el()
        idlist = mesh.slice(plane,component=[(2,1)])
        np.testing.assert_equal(
                                mesh.comp.values[mesh.xyz_el.values[:,0] <.1],
                                1
                                )
        np.testing.assert_equal(
                                mesh.comp.values[mesh.xyz_el.values[:,0] >.1],
                                2
                                )

            #steam.io.write_vtu(name="test.vtu",mesh=mesh)

    def test_returnids( self ):
        """ Slice and check the returned ids."""
        plane = np.zeros((1,3,3))
        plane[0,1,1] = 1
        plane[0,2,2] = 1
        
        mesh = steam.mesh.mk_cube( 1, 5)
        mesh.get_xyz_el()
        o_mesh = steam.mesh.Mesh(copy=mesh)
        o_soln = steam.solution.uniform_soln(mesh,values=['rand'])
        
        idlist = mesh.slice(plane,returnList=True)
        soln = steam.solution.uniform_soln(mesh)
        # Use the idlist to expand the solution:
        soln.data['q1'] = o_soln.data.values[idlist]
        
        # Since we didn't cut, this shoudl be one-to-one
        np.testing.assert_equal(
                                soln.data.values,
                                o_soln.data.values
                                )
        
        ### Now do some cutting
        plane[0,:,0] = .1
        mesh = steam.mesh.mk_cube( 1, 5)
        mesh.get_xyz_el()
        
        idlist = mesh.slice(plane,returnList=True)
        soln = steam.solution.uniform_soln(mesh)
        # Use the idlist to expand the solution:
        soln.data['q1'] = o_soln.data.values[idlist]
        
        # Now, let's compare to a closest-point transformation
        trans = steam.interpolate.Transform(o_mesh,mesh)
            # I have to do them all individually so that we don't snap to the
            # closest point in an adjacent face
        trans.inverse_distance(k=1,source_comp=['1'],target_comp=['1'])
        trans.inverse_distance(k=1,source_comp=['2'],target_comp=['2'])
        trans.inverse_distance(k=1,source_comp=['3'],target_comp=['3'])
        trans.inverse_distance(k=1,source_comp=['4'],target_comp=['4'])
        trans.inverse_distance(k=1,source_comp=['5'],target_comp=['5'])
        trans.inverse_distance(k=1,source_comp=['6'],target_comp=['6'])
        check_soln = trans.apply(o_soln)
        
        # Since we didn't cut, this shoudl be one-to-one
        np.testing.assert_equal(
                                check_soln.data.values,
                                soln.data.values
                                )

        #steam.io.write_vtu(name="test.vtu",mesh=mesh)

    @classmethod
    def tearDownClass( cls ):
        """Clean out temporary files."""

        import os

        files_to_remove = ( 'test.vtu' )
        for f in files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )


pass
        
        
#! Run the tests
if __name__ == '__main__':
    ut.main()

