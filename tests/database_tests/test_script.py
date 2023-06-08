#! test_script.py
#! In order to run these tests with verbose output:
#! >>> python test_script.py -v
#!
"""
Simple script that makes a database in 5D space,
creates a database object and loads map varaibles,
and then interpolates to 50 points inside the space.

Linear variables are compared to truth values on the
maps and dependant variables in the database.  A
difference between the truth and interpolated values
results in an error code.

Database reading and writing are tested as are
the accuracy of the interpolation results.

Note: This is just a test of point solution data.
"""

import scipy as sp
import numpy as np
import unittest as ut
import pandas as pd

import os
import steam.database  as dbms_db
import steam.table     as dbms_tab
import steam.container as dbms_cont

class TestDatabase( ut.TestCase ):


    #! 5-D interpolation inputs.  Can up ndim to as high as 10 (>8 takes forever).
    ndim = 5
    coef = np.array([.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0])
    vars =  ['Mach','alpha','beta','rho','temp','foo','bar','test','ten','dims']

    #! This adjusts the size of the meshes
    npts  = 100     # ~6MB (with 1024 fs_pts)
    #npts  = 1000   # 62MB (with 1024 fs_pts)
    #npts  = 20000  #1.6GB (with 1024 fs_pts)

    fs_pts = 1024
    fs_pts = 50
    
    #! Going to need 2 points per dimension, so 2**ndim points
    bounds = np.zeros((ndim,2))
    bounds[:,0] = 0.0
    bounds[:,1] = 1.0

    #! Fix the random number seed
    np.random.seed(0)

    #!< Create the necessary sample files that will be used to test 
    #!  interpolation
    #!>
    @classmethod
    def setUpClass( self ):
        """This method gets called before any of the tests are run in 
        order to create the fake data that is needed.  Large files can
        be made by adjusting the mesh point counts (above).

        This goes ahead and makes the HDF5 file.
        """ 

        from itertools import product  # This is super cool

        #! First make the database points by using the bounds
        points = np.array([]).reshape(0,self.ndim)
        for set in product(*self.bounds):
            points = np.vstack((points,set))
        
        #! Add interior Points (Randomly!)
        for i in range(self.fs_pts):
            point = np.random.rand(1,self.ndim)
            points = np.vstack((points,point))
        
        #! Add values
        values = np.array([]).reshape(0,2)
        for point in points:
            vect = point                            # Linear
            line = np.dot(vect,self.coef[0:self.ndim])
            vect = point * point                    # Quadratic
            quad = np.dot(vect,self.coef[0:self.ndim])
            # Have to add this as a tuple.
            values = np.vstack((values,[line,quad]))

        #! Create the mesh as a 2D region in 3D space
        mesh      = np.random.rand(self.npts,3) *2.0-1.0  
        mesh[:,2] = 0.0
        gdf = pd.DataFrame(mesh, columns = ['X','Y','Z'])
        gdf.to_csv('mesh.cdat',index=False,sep=' ')
        
        #! A solution is going to have a simple equation with four solutions at each point:
        #! - q1 = coef*point   - sqrt(X^2+Y^2)
        #! - q2 = coef*point   -     (X^2+Y^2)
        #! - q3 = coef*point^2 - sqrt(X^2+Y^2)
        #! - q4 = coef*point^2 -     (X^2+Y^2)
        
        #! Make the solution directory
        if not os.path.isdir('sols'):
            os.mkdir('sols')

        sols = np.array([]).reshape(0,2)
        for (i,point) in enumerate(points):
            vect = point                            # Linear
            line = np.dot(vect,self.coef[0:self.ndim])
            vect = point * point                    # Quadratic
            quad = np.dot(vect,self.coef[0:self.ndim])
        
            path = "sols/{:06d}.cdat".format(i)
        
            # Have to add this as a tuple.  type and path
            sols = np.vstack((sols,["DF",path]))
        
            #! Make the solution
            q = np.array([]).reshape(0,4)
            for xyz in mesh:
                q1 = line - np.dot(xyz,xyz)**0.5
                q2 = line - np.dot(xyz,xyz)
                q3 = quad - np.dot(xyz,xyz)**0.5
                q4 = quad - np.dot(xyz,xyz)
                q = np.vstack((q,[q1,q2,q3,q4]))
        
            qdf = pd.DataFrame(q, columns = ['q1','q2','q3','q4'])
        
            qdf.to_csv(path,index=False,sep=' ')
        
        #! Make this a dataframe and write it out
        df = np.hstack((points,values,sols))
        dfvars = (self.vars[0:self.ndim]).copy()
        dfvars.extend(['linear','quad'])
        dfvars.extend(['soln_type' ,'soln_path'])
        
        dff = pd.DataFrame(df, columns = dfvars)
        dff.to_csv('db.cdat',index=False,sep=' ')

        #! Make the database from the sample data.
        test_db = dbms_db.Database()
        test_db.read_cdat('db.cdat')
        test_db.read_mesh('mesh.cdat','DF')
        test_db.load_soln()

        #! Make the fake database
        target_df = np.array([]).reshape(0,self.ndim)

        #! Add target points (Randomly!)
        for i in range(50):
            point     = np.random.rand(1,self.ndim)
            target_df = np.vstack((target_df,point))

        
        target = dbms_tab.Table()
        target.data = pd.DataFrame(target_df, columns = self.vars[0:self.ndim])

        cont = dbms_cont.Container("test.h5")
        cont.add("/database",test_db)
        cont.add("/target" ,target)
        cont.write()

 
    def test_iter( self ):

        #! Load and unpack the container
        cont = dbms_cont.Container("test.h5")
        cont.read()
        db = cont.obj['database']

        for (i,p) in enumerate(db):
            pass

        self.assertEqual(i,len(db.data)-1)

        p2 = db.get_point(len(db.data)-1)
        self.assertTrue(np.all(p.data == p2.data))
        self.assertTrue(p.mesh == p2.mesh)
        self.assertTrue(p.soln == p2.soln)

        for p in db:
            pass

        self.assertTrue(np.all(p.data == p2.data))
        self.assertTrue(p.mesh == p2.mesh)
        self.assertTrue(p.soln == p2.soln)
 
    def test_get_add( self ):

        #! Load and unpack the container
        cont = dbms_cont.Container("test.h5")
        cont.read()
        database = cont.obj['database']

        loc = 5
        point = database.get_point(loc)

        glen = len(database.meshes)
        slen = len(database.solns)

        for key in (point.data.keys()):
            self.assertEqual(
                        point.data[key],
                        database.data[key][loc]
                        )

        ### This should not add a new mesh since
        ### it is a copy of an existing one.
        ### Same for solution.
        database.add_point(point)
        self.assertEqual(
                    len(database.meshes),
                    glen
                    )

        ### This should not add a new soln since
        ### it is a copy of an existing one
        self.assertEqual(
                    len(database.solns),
                    slen
                    )

    #!< Make a target database and interpolate to it.
    #!  Since this is linear interpolation, we only
    #!  expect the linear variable to agree, the
    #!  quad variable is just along for the ride.
    #!>
    def test_interpolation( self ):

        #! Load and unpack the container
        cont = dbms_cont.Container("test.h5")
        cont.read()
        database = cont.obj['database']
        target   = cont.obj['target']

        database.set_indep(self.vars[0:self.ndim])
        database.set_dep(['linear','quad'])

        database.new_intspace()

        #! Collect target cases into dataframe
        cases    = target.data[self.vars[0:self.ndim]]

        #! Do the interpolation
        DBPoints = database.interpolate(cases)

        #! Interpolate to data and compare with exact results
        for (i,db_point) in enumerate(DBPoints):

            case = cases.iloc[i].values

            #! Compare the data table results
            vect = case
            line = np.dot(case,self.coef[0:self.ndim])
            vect = case * case                     # Quadratic
            quad = np.dot(vect,self.coef[0:self.ndim])

            self.assertEqual(
                            "{:12e}".format(db_point.data['linear']),
                            "{:12e}".format(line)
                            )

            #! Old Items that would print out the error
            #error = abs(data['linear']-line)
            #print("   Error     : {:.12f}".format(error))
            #print("   Quad Diff : {:.12f}".format(data['quad']  -quad))

            #! Now run through each solution
            soln_error = {'q1':0.0,'q2':0.0,'q3':0.0,'q4':0.0}
            
            #! A solution is going to have a simple equation with four solutions at each point:
            #! - q1 = coef*point   - sqrt(X^2+Y^2)
            #! - q2 = coef*point   -     (X^2+Y^2)
            #! - q3 = coef*point^2 - sqrt(X^2+Y^2)
            #! - q4 = coef*point^2 -     (X^2+Y^2)

            for (i,row) in db_point.soln.data.iterrows():
                xyz = db_point.mesh.xyz_pt.iloc[i]
                q1 = line - np.dot(xyz,xyz)**0.5
                q2 = line - np.dot(xyz,xyz)
                q3 = quad - np.dot(xyz,xyz)**0.5
                q4 = quad - np.dot(xyz,xyz)
                soln_error['q1'] += abs(row['q1'] - q1)
                soln_error['q2'] += abs(row['q2'] - q2)
                soln_error['q3'] += abs(row['q3'] - q3)
                soln_error['q4'] += abs(row['q4'] - q4)

                self.assertEqual(
                                 "{:12e}".format(row['q1']) ,
                                 "{:12e}".format(     q1)
                                )
                self.assertEqual(
                                 "{:12e}".format(row['q2']) ,
                                 "{:12e}".format(     q2)
                                )

            #! Old items that would print out the error
            #print("   q1 Error  : {:.12f}".format( soln_error['q1']))
            #print("   q2 Error  : {:.12f}".format( soln_error['q2']))
            #print("   q3 Diff   : {:.12f}".format( soln_error['q3']))
            #print("   q4 Diff   : {:.12f}".format( soln_error['q4']))

            #stol = 1.0e-15
            #mtol = 1.0e-15*len(soln)

            #if (error > stol or soln_error['q1'] > mtol or soln_error['q2'] > mtol):
            #    print(stol,mtol,mtol)
            #    print(error,soln_error['q1'],soln_error['q2'])
            #    print(" !! ** Error in interpolation ** !!")
            #sys.exit(1)


    def test_remove_pt( self ):
        """ Remove boundary points from database and interpolate.

        Confirm that the results don't change in the interior.

        One of the main challenges of this test, and one that is somewhat subtle,
        is that sub_db1 and sub_db2 have the same data stored in <db>.data, except
        that the indices differ.  This is the differences that exercises the 
        interpolation code and its change from the use of DataFrame.iloc to 
        DataFrame.loc.
        """

        ### Load and unpack the container
        with dbms_cont.Container("test.h5", 'r' ) as cont:
            cont.read()
            db = cont.obj['database']
            target   = cont.obj['target']

        db.set_indep(self.vars[0:self.ndim])
        db.set_dep(['linear','quad'])

        ### Create a database with a subset of the points in db
        #   Populate it by re-reading db and removing points with Mach <= 0.5
        with dbms_cont.Container("test.h5", 'r' ) as cont:
            cont.read()
            sub_db1 = cont.obj['database']
        for pt in sub_db1:
            if pt.data.Mach <= 0.5:
                sub_db1.remove_pt( pt )
        sub_db1.set_indep(self.vars[0:self.ndim])
        sub_db1.set_dep(['linear','quad'])

        ### Create another database with a subset of the points in db
        #   Populate it by adding DBPoints from db if Mach > 0.5
        sub_db2 = dbms_db.Database()
        for pt in db:
            if pt.data.Mach > 0.5:
                sub_db2.add_point( pt )
        sub_db2.set_indep(self.vars[0:self.ndim])
        sub_db2.set_dep(['linear','quad'])

        ### Define interpolation spaces for all databases
        db.new_intspace('high_m', subset='Mach > 0.5')
        sub_db1.new_intspace( 'entire' )
        sub_db2.new_intspace( 'entire' )

        ### Collect target cases into dataframe
        #   Drop cases with Mach <= 0.5
        cases    = target.data[self.vars[0:self.ndim]].query( 'Mach > 0.5' )

        ### Perform interpolation in all 3 databases
        pts0 = db.interpolate(     cases, 'high_m', quiet=True, return_empty=False)
        pts1 = sub_db1.interpolate(cases, 'entire', quiet=True, return_empty=False)
        pts2 = sub_db2.interpolate(cases, 'entire', quiet=True, return_empty=False)

        self.assertEqual( len(pts0), len(pts1) )
        self.assertEqual( len(pts1), len(pts2) )

        for pt0, pt1, pt2 in zip( pts0, pts1, pts2 ):
            ### Test the data
            self.assertTrue( pt0.data.equals(pt1.data) )
            self.assertTrue( pt1.data.equals(pt2.data) )

            ### Test the meshes
            self.assertEqual( pt0.mesh, pt1.mesh )
            self.assertEqual( pt1.mesh, pt2.mesh )

            ### Test the solutions
            self.assertEqual( pt0.soln, pt1.soln )
            self.assertEqual( pt1.soln, pt2.soln )

    #!< Just remove all of the files that we made in the
    #! setUpClass routine.
    #!>
    @classmethod
    def tearDownClass( self):
        """Clean out temporary files made as well as the
        HDF5 file.
        """
        
        os.remove('db.cdat')
        os.remove('mesh.cdat')
        os.remove('test.h5')
        import shutil
        shutil.rmtree('sols')

#! Run the tests
if __name__ == '__main__':
    ut.main()
