#! /bin/env python

import unittest as ut
import numpy as np
import pandas as pd
import steam

def get_simple_cube(np,q1):
    mesh = steam.mesh.mk_cube( 1, np)
    mesh.get_xyz_el()
    soln = steam.solution.uniform_soln( mesh, values = [q1])
    return (mesh,soln)
# Does it depend on LibMesh?
class TestProject( ut.TestCase ):

    def init(self):
        # Ok, compare some grids.  Make a database of 3 grids and then
        # replace them with one of their own.
        db = steam.database.Database()
        dbp = steam.database.DBPoint()

        for np in [5,6,7,6,5]:
            dbp.data = pd.Series({'np':np})
            (dbp.mesh,dbp.soln) = get_simple_cube(np,np)
            db.add_point(dbp)
                    
        cont = steam.container.Container("mem.h5",'w')
        cont.add("/database",db)
        cont.write()
        cont = steam.container.Container("disk.h5",'w')
        cont.add("/database",db)
        cont.write()

    def do_bothid_tests(self,db):
        """ Simple routine that will test update_meshid and update_solnid.

        Important that the input has three different meshes/solns and
        that the output collapses it down to one of each.  Have fun
        in between.
        """

        # These should just hang out
        for i in range(len(db.data)):
            point  = db.get_point(i)
            old_id = point.get_meshid()
            new_id = db.update_meshid(old_id,point.mesh)
            self.assertTrue(old_id == new_id)
        self.assertTrue(len(db.meshes) == 3)

        # Update all of the solutions to use the last mesh
        for gkey in list(db.meshes.keys()):
            db.update_meshid(gkey,point.mesh)
         
        self.assertTrue(len(db.meshes) == 1)
        self.assertTrue(len(db.solns)  == 3)

        # Replace them all at once:
        (mesh,soln) = get_simple_cube(3,3)
        db.update_meshid(list(db.meshes.keys())[0],mesh)
        self.assertTrue(len(db.meshes) == 1)


        print(" - Now do solutions")
        # Now try and update the solution with redundant data
        for i in range(len(db.data)):
            point  = db.get_point(i)
            old_id = point.get_solnid()
            new_id = db.update_solnid(old_id,point.soln)
            self.assertTrue(old_id == new_id)
        self.assertTrue(len(db.solns) == 3)

        # Update all of the solutions to use the last mesh
        for sid in list(db.solns.keys()):
            db.update_solnid(sid,point.soln)
        self.assertTrue(len(db.solns)  == 1)

        # Replace them all at once:
        (mesh,soln) = get_simple_cube(3,3)
        db.update_solnid(list(db.solns.keys())[0],soln)
        self.assertTrue(len(db.solns) == 1)

        print(" - Raise all of the exceptions")
        mid = list(db.meshes.keys())[0]
        sid = list(db.solns.keys())[0]
        # Now, raise all of the exceptons
        # - Not the correct type
        with self.assertRaisesRegex(TypeError,"not a STEAM"):
            db.update_meshid(mid,"string")
        with self.assertRaisesRegex(TypeError,"not a STEAM"):
            db.update_solnid(sid,"string")
        # - Not a correct id
        with self.assertRaisesRegex(LookupError, ".*not found in database.*"):
            db.update_meshid("foobar",mesh)
        with self.assertRaisesRegex(LookupError, ".*not found in database.*"):
            db.update_solnid("foobar",soln)

    def test_update_ids( self ):
        """ Test the database.update_solnid and database.update_meshid methods."""

        self.init()

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        disk_cont = steam.container.Container('disk.h5','a')
        disk_cont.read(options={'database':{'onDisk':True}})
        disk_db   = disk_cont.obj['database']

        self.do_bothid_tests(mem_db)
        self.do_bothid_tests(disk_db)

        # Write the onDisk database out and read it in.  Make
        # sure that it is the correct size when read back in.
        disk_db.onDisk_to_mem()
        cont = steam.container.Container("disk_post1.h5",'w')
        cont.add("database",disk_db)
        cont.write()

        disk_p1_cont = steam.container.Container('disk_post1.h5')
        disk_p1_cont.read()
        disk_p1_db   = disk_p1_cont.obj['database']
        self.assertTrue(len(disk_p1_db.meshes) == 1)
        self.assertTrue(len(disk_p1_db.solns)  == 1)


    def do_bothidx_tests(self,db):
        """ Simple routine that will test update_meshidx and update_solnidx.

        Important that the input has three different meshes/solns and
        that the output collapses it down to one of each.  Have fun
        in between.
        """

        # These should just hang out
        for i in range(len(db.data)):
            point  = db.get_point(i)
            old_id = point.get_meshid()
            new_id = db.update_meshidx(i,point.mesh)
            self.assertTrue(old_id == new_id)
        self.assertTrue(len(db.meshes) == 3)

        # Update all of the solutions to use the last mesh
        for i in range(len(db.data)):
            db.update_meshidx(i,point.mesh)
         
        self.assertTrue(len(db.meshes) == 1)
        self.assertTrue(len(db.solns)  == 3)


        print(" - Now do solutions")
        # Now try and update the solution with redundant data
        for i in range(len(db.data)):
            point  = db.get_point(i)
            old_id = point.get_solnid()
            new_id = db.update_solnidx(i,point.soln)
            self.assertTrue(old_id == new_id)
        self.assertTrue(len(db.solns) == 3)

        # Update all of the solutions to use the last mesh
        for i in range(len(db.data)):
            db.update_solnidx(i,point.soln)
        self.assertTrue(len(db.solns)  == 1)

        (mesh,soln) = get_simple_cube(3,3)

        print(" - Raise all of the exceptions")
        mid = list(db.meshes.keys())[0]
        sid = list(db.solns.keys())[0]
        # Now, raise all of the exceptons
        # - Not the correct type
        with self.assertRaisesRegex(TypeError,"not a STEAM"):
            db.update_meshid(mid,"string")
        with self.assertRaisesRegex(TypeError,"not a STEAM"):
            db.update_solnid(sid,"string")
        # - Not a correct idx
        with self.assertRaisesRegex(LookupError, ".*not find .* at index.*"):
            db.update_meshidx(10,mesh)
        with self.assertRaisesRegex(LookupError, ".*not find .* at index.*"):
            db.update_solnidx(20,soln)




    def test_update_idxs( self ):
        """ Test the database.update_solnidx and database.update_meshidx methods.
        """

        self.init()

        mem_cont  = steam.container.Container('mem.h5')
        mem_cont.read()
        mem_db    = mem_cont.obj['database']

        disk_cont = steam.container.Container('disk.h5','a')
        disk_cont.read(options={'database':{'onDisk':True}})
        disk_db   = disk_cont.obj['database']

        self.do_bothidx_tests(mem_db)
        self.do_bothidx_tests(disk_db)

        # Write the onDisk database out and read it in.  Make
        # sure that it is the correct size when read back in.
        disk_db.onDisk_to_mem()
        cont = steam.container.Container("disk_post1.h5",'w')
        cont.add("database",disk_db)
        cont.write()

        disk_p1_cont = steam.container.Container('disk_post1.h5')
        disk_p1_cont.read()
        disk_p1_db   = disk_p1_cont.obj['database']
        self.assertTrue(len(disk_p1_db.meshes) == 1)
        self.assertTrue(len(disk_p1_db.solns)  == 1)

    @classmethod
    def tearDownClass( cls ):
        """This method gets called once, at the very end of testing in order
        to clean out temporary files.
        """
        import os
 
        files_to_remove = ['disk.h5','mem.h5','disk_post1.h5']
        for f in files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )


#! Run the tests
if __name__ == '__main__':
    ut.main()
