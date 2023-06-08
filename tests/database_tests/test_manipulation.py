#! /bin/env python

import unittest as ut
import numpy as np
from numpy.random import rand
import pandas as pd
import steam
from code import interact

def get_simple_cube(q1):
    mesh = steam.mesh.mk_cube( 1, 5)
    mesh.get_xyz_el()
    soln = steam.solution.uniform_soln( mesh, values = [q1])
    return (mesh,soln)

# Does it depend on LibMesh?
class TestProject( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """Create a 2D database of cubes."""

        db = steam.database.Database()
        dbp = steam.database.DBPoint()

        # Database is linear, just the product of v1 and v2
        print("START")
        import time
        for v1 in range(5):
            start = time.time()
            for v2 in range(5):
                dbp.data = pd.Series({'v1':v1,'v2':v2})
                (dbp.mesh,dbp.soln) = get_simple_cube(v1*v2)
                db.add_point(dbp)
            print(v1,len(db.data),time.time() - start)
            print("   ",len(db.meshes),len(db.solns))
                
        cont = steam.container.Container("mem.h5",'w')
        cont.add("/database",db)
        cont.write()
        cont = steam.container.Container("disk.h5",'w')
        cont.add("/database",db)
        cont.write()

    def test_eq( self ):
        """ Test functionality of the equality operator, "=="

        As of 4/29/2020, no proper __eq__ function has been written.  
        All that is tested, therefore, is that the proper error is raised.
        """
    
        db1 = steam.database.Database()
        db2 = steam.database.Database()
        with self.assertRaises( NotImplementedError ) as cm:
            db1 == db1

        with self.assertRaises( NotImplementedError ) as cm:
            db1 == db2

        with self.assertRaises( NotImplementedError ) as cm:
            db2 == db1
            
        with self.assertRaises( NotImplementedError ) as cm:
            db1.__eq__( db2 )

    def test_formats( self ):
        """Add a new variable to the database with different input formats"""

        ### Load the database
        with steam.container.Container( 'mem.h5', 'r' ) as cont:
            cont.read()
        db = cont.obj['database']

        ### New data
        new_dat = rand( db.data.shape[0] )

        ### Add data in multiple formats
        db.add_var( 'list',   list(new_dat) )
        db.add_var( 'tuple',  tuple(new_dat) )
        db.add_var( 'array',  np.array(new_dat) )
        db.add_var( 'series', pd.Series(new_dat) )

        ### Check that all added columns are equal
        self.assertTrue( db.data.list.equals( 
                         db.data.tuple) )
        self.assertTrue( db.data.tuple.equals( 
                         db.data.array) )
        self.assertTrue( db.data.array.equals( 
                         db.data.series) )
        self.assertTrue( db.data.series.equals( 
                         db.data.list) )

    def test_add_with_index( self ):
        """ Use add_var with optional index_in input """

        ### Load the database
        with steam.container.Container( 'mem.h5', 'r' ) as cont:
            cont.read()
        db = cont.obj['database']

        ### New data
        new_dat = rand( db.data.shape[0] )
        
        ### Add data with the same index as the solution
        db.add_var( 'full_dat', new_dat, 
                                index_in=db.data.index )
        for val in db.data.full_dat:
            self.assertFalse( np.isnan( val ) )
        
        ### Add discontinuous data and make sure it appears in the right place
        disc_dat = list( range(0, db.data.shape[0], 2) )
        db.add_var( 'discontinuous', disc_dat, index_in=disc_dat )
        for ind, val in db.data.discontinuous.items():
            if ind in disc_dat:
                self.assertEqual( val, ind )
            else:
                self.assertTrue( np.isnan( val ) )

    def test_get_soln( self ):
        """Make sure that all ways to get_soln return a copy that can be manipulated 
        without changing the underlying data.
        """

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        disk_cont = steam.container.Container('disk.h5','a')
        disk_cont.read(options={'database':{'onDisk':True}})
        disk_db   = disk_cont.obj['database']

        #! Memory
        value = mem_db.data['v1'][0] * mem_db.data['v2'][0]
        sid   = mem_db.data['solnid'][0]

        copy = mem_db.get_soln(sid)
        copy.data['q1'] = 999.0
            # Did it change as intended
        np.testing.assert_allclose(copy.data.values,999.0)
            # One more test
        (mesh,test) = get_simple_cube(999)
        np.testing.assert_allclose(copy.data.values,test.data.values)
            # The do no harm test - is the database copy the same as it was?
        (mesh,soln) = get_simple_cube(value)
        np.testing.assert_allclose(mem_db.solns[sid].data.values,soln.data.values)
 
        #! Disk
        value = mem_db.data['v1'][1] * mem_db.data['v2'][1]
        sid   = mem_db.data['solnid'][1]

        copy = disk_db.get_soln(sid)
        copy.data['q1'] = 999.0
            # Did it change as intended
        np.testing.assert_allclose(copy.data.values,999.0)
            # One more test
        (mesh,test) = get_simple_cube(999)
        np.testing.assert_allclose(copy.data.values,test.data.values)
            # The do no harm test - is the database copy the same as it was?
        (mesh,soln) = get_simple_cube(value)
        disk_soln = disk_db.get_soln(sid)
        np.testing.assert_allclose(disk_soln.data.values,soln.data.values)

    def test_get_mesh( self ):
        """Make sure that get_mesh handles copies correctly.
        """

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        disk_cont = steam.container.Container('disk.h5','a')
        disk_cont.read(options={'database':{'onDisk':True}})
        disk_db   = disk_cont.obj['database']


        mid   = mem_db.data['meshid'][0]

        #! Memory
        copy = mem_db.get_mesh(mid)
        copy.transform(["t10x","t10y","t10z"])
            # Did it change as intended
        np.testing.assert_allclose(copy.xyz_pt,
                                   mem_db.meshes[mid].xyz_pt+10.0)
            # The do no harm test - is the database copy the same as it was?
        (mesh,soln) = get_simple_cube(0.0)
        self.assertTrue(mesh == mem_db.meshes[mid])
 

        #! Disk 
        copy = disk_db.get_mesh(mid)
        copy.transform(["t10x","t10y","t10z"])
            # Did it change as intended
        np.testing.assert_allclose(copy.xyz_pt,
                                   disk_db.meshes[mid].xyz_pt+10.0)
            # The do no harm test - is the database copy the same as it was?
        (mesh,soln) = get_simple_cube(0.0)
        self.assertTrue(mesh == disk_db.meshes[mid])

    def test_update_meshid( self ):
        """Make sure that update_meshid works.
        """

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        disk_cont = steam.container.Container('disk.h5','a')
        disk_cont.read(options={'database':{'onDisk':True}})
        disk_db   = disk_cont.obj['database']


        mid   = mem_db.data['meshid'][0]

        #! Memory
        copy = mem_db.get_mesh(mid)
        copy.transform(["t10x","t10y","t10z"])
            # Did it change as intended
        np.testing.assert_allclose(copy.xyz_pt,
                                   mem_db.meshes[mid].xyz_pt+10.0)
        # Update the database
        mem_db.update_meshid(mid,copy)
        self.assertTrue(copy == mem_db.meshes[mid])
 

        #! Disk 
        copy = disk_db.get_mesh(mid)
        copy.transform(["t10x","t10y","t10z"])
            # Did it change as intended
        np.testing.assert_allclose(copy.xyz_pt,
                                   disk_db.meshes[mid].xyz_pt+10.0)
        
        # Update the database
        disk_db.update_meshid(mid,copy)
        self.assertTrue(copy == disk_db.meshes[mid])

    

    def test_update_solnid( self ):
        """Make sure that update_solnid works.
        """

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        disk_cont = steam.container.Container('disk.h5','a')
        disk_cont.read(options={'database':{'onDisk':True}})
        disk_db   = disk_cont.obj['database']

        #! Memory 
        value = mem_db.data['v1'][2] * mem_db.data['v2'][2]
        sid   = mem_db.data['solnid'][2]

        copy = mem_db.get_soln(sid)
        copy.data['q1'] = 999.0
            # Did it change as intended
        np.testing.assert_allclose(copy.data.values,999.0)
            # One more test
        (mesh,test) = get_simple_cube(999)
        np.testing.assert_allclose(copy.data.values,test.data.values)

        # Update the database
        mem_db.update_solnid(sid,copy)
            # This time, the database should have changed
        np.testing.assert_allclose(mem_db.solns[sid].data.values,copy.data.values)
 

        #! Disk
        value = mem_db.data['v1'][3] * mem_db.data['v2'][3]
        sid   = mem_db.data['solnid'][3]

        copy = disk_db.get_soln(sid)
        copy.data['q1'] = 999.0
            # Did it change as intended
        np.testing.assert_allclose(copy.data.values,999.0)
            # One more test
        (mesh,test) = get_simple_cube(999)
        np.testing.assert_allclose(copy.data.values,test.data.values)

        # Update database
        disk_db.update_solnid(sid,copy)
            # This time, the database should have changed
        disk_soln = disk_db.get_soln(sid)
        np.testing.assert_allclose(disk_soln.data.values,
                                        copy.data.values)

    def test_update_point( self):
        """Make sure that update_mesh and update_soln work for a DBPoint.
        """

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        disk_cont = steam.container.Container('disk.h5','a')
        disk_cont.read(options={'database':{'onDisk':True}})
        disk_db   = disk_cont.obj['database']

        #! Memory
        point = mem_db.get_point(7)
        mid = point.get_meshid()
        sid = point.get_solnid()

        point.mesh.transform(["t10x","t10y","t10z"])
        point.soln.data['q1'] = 999.0

        # Did it change as intended
        np.testing.assert_allclose(point.mesh.xyz_pt,
                                   mem_db.meshes[mid].xyz_pt+10.0)
        np.testing.assert_allclose(point.soln.data.values,999.0)

        # Update the database
        point.update_mesh()
        self.assertTrue(point.mesh == mem_db.meshes[point.get_meshid()])

        point.update_soln()
        np.testing.assert_allclose(mem_db.solns[point.get_solnid()].data.values,
                                   point.soln.data.values)
 

        #! Disk
        point = disk_db.get_point(8)
        mid = point.get_meshid()
        sid = point.get_solnid()

        point.mesh.transform(["t10x","t10y","t10z"])
        point.soln.data['q1'] = 999.0

        # Did it change as intended
        np.testing.assert_allclose(point.mesh.xyz_pt,
                                   disk_db.meshes[mid].xyz_pt+10.0)
        np.testing.assert_allclose(point.soln.data.values,999.0)

        # Update the database
        point.update_mesh()
        self.assertTrue(point.mesh == disk_db.meshes[point.get_meshid()])

        point.update_soln()
        soln = disk_db.get_soln(point.get_solnid())
        np.testing.assert_allclose(      soln.data.values,
                                   point.soln.data.values)


    def test_update_point_loop( self):
        """Test update_mesh and update_soln for a database point in a loop..
        """

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        #! Memory
        num_before = len(list(set(mem_db.data.meshid)))
        for (i,point) in enumerate(mem_db):
            point.mesh.transform(["t{}x".format(i)])
            point.update_mesh()
        num_after  = len(list(set(mem_db.data.meshid)))

        self.assertTrue(num_before == 1)
        self.assertTrue(num_after  == mem_db.data.shape[0])

        for (i,point) in enumerate(mem_db):
            # Since the original data is a linear function
            # that loosely corresponds to the index, muck
            # things up to ensure we don't accidentally match
            point.soln.data += 1.5
            point.soln.data *= (i+2)**2.0
            point.update_soln()
        num_after  = len(list(set(mem_db.data.solnid)))
        self.assertTrue(num_after  == mem_db.data.shape[0])

    def test_loop_add( self ):
        """Confirm that adding points in a loop doesn't go on infinitely."""

        with steam.container.Container( 'mem.h5', 'r' ) as cont:
            cont.read()
            db = cont.obj['database']

        original_len = len( db.data )

        for pt_count, pt in enumerate( db ):
            if pt_count >= original_len:
                self.fail( 'Getting into infinite add loop.' )
            pt.data.v1 += 4*original_len
            pt.data.v2 += 6*original_len
            db.add_point( pt )

        self.assertEqual( len(db.data), 2*original_len )

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """
        import os

        files_to_remove = ['test.h5','disk.h5','mem.h5']
        for f in files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )


#! Run the tests
if __name__ == '__main__':
    ut.main()
