#! /bin/env python

import unittest as ut
import os
import steam
from copy import deepcopy
from code import interact

class TestContainer( ut.TestCase ):

    @classmethod
    def setUpClass( cls ):
        """ Create the list of files to delete """
        cls.files_to_remove = []

    def setUp( self ):
        """This method gets called before each test_* method in order to 
        prepare the proper variables/data.
        """ 
        # I just want to create some stock thing that we can add to
        # containers:
        self.mesh   = steam.mesh.mk_cube(1.0,5)
        self.tab    = steam.table.Table('aero.cdat')
        self.data   = steam.database.Database('aero.cdat')
        self.data.add_mesh(self.mesh)
        self.string = "This is a test string, and it is very important to me."

        # List of objects with bad characters in their names
        # Bad means non-HDF5 conforming: does not match the pattern ``^[a-zA-Z_][a-zA-Z0-9_]*$``
        self.bad_chars = ["bad/1_test",
                          "bad/1.0test",
                          "bad/test_1.0",
                          "bad/this-that",
                          "bad/this(that)",
                          "bad/_foo_",
                          "bad/%this"]
        pass

    def tearDown(self):
        """ Just remove the test HDF5 file"""
        if os.path.isfile("test1.h5"):
            os.remove("test1.h5")

        for f in self.files_to_remove:
            if os.path.isfile( f ):
                os.remove( f )

    def test_read_and_write( self ):
        """ Simple test to read/write to/from a container
        """

        self.files_to_remove.append( 'test1.h5' )

        name = "Test Container"
        with steam.container.Container("test1.h5",
                                       desc=name
                                       ) as c:
            c.add("mesh",self.mesh)
            c.add("table",self.tab)
            c.add("data",self.data)
            c.add("str",self.string)
            for var in self.bad_chars:
                c.add(var,var)

        # Now try and read from it
        with steam.container.Container("test1.h5") as c:
            c.read()
            mesh2   = c.obj["mesh"]
            table2  = c.obj["table"]
            data2   = c.obj["data"]
            string2 = c.obj["str"]
            name2   = c.description
            steam_version = c.obj['steam_version']

        self.assertTrue(name2   == name)
        self.assertTrue(mesh2   == self.mesh)
        self.assertTrue(data2.data.equals(self.data.data))
        self.assertTrue(string2 == self.string)
        self.assertTrue(table2  == self.tab)
        self.assertEqual( steam_version, steam.__version__ )

        # Check the bad vars
        for var in self.bad_chars:
            self.assertTrue(c.obj[var]   == var)

        pass

    def test_old_version( self ):
        """ Fake that the container was written with an older version of STEAM
        """

        name = "Test Container"
        old_container = 'test_old.h5'
        self.files_to_remove.append( old_container )
        with steam.container.Container( old_container,
                                       desc=name
                                       ) as c:
            c.add("mesh",self.mesh)
            c.add( 'steam_version', 'not_a_version' )

        ### Read the "old container"
        with steam.container.Container( old_container ) as c2:
            c2.read()

        self.assertNotEqual( steam.__version__, c2.obj['steam_version'] )

    def test_change_type(self):
        """ Test what happens when you change the type of an object
        between adding it and writing. """

        name = "Test Container"
        with steam.container.Container("test1.h5",
                                       desc=name
                                       ) as c:
            c.add("mesh",self.mesh)
            c.obj["mesh"] = self.data

        # Now try and read from it
        with steam.container.Container("test1.h5") as c:
            c.read()
            mesh2   = c.obj["mesh"]

        self.assertFalse(mesh2   == self.mesh)
        self.assertTrue(mesh2.data.equals(self.data.data))

    def test_write_subset(self):
        """ Just write a subset of the objects and make sure they check out."""


        name = "Test Container"
        cout = steam.container.Container("test1.h5", desc=name)
        cout.add("mesh",self.mesh)
        cout.add("table",self.tab)
        cout.add("data",self.data)
        cout.add("str",self.string)
        cout.write(items = ["mesh","str"])

        cin  = steam.container.Container("test1.h5")
        cin.read()
        self.assertTrue(cin.obj["mesh"] == self.mesh)
        self.assertTrue(cin.obj["str" ] == self.string)

    def test_delete_subset(self):
        """ Delete a subset of a container."""

        name = "Test Container"
        cout = steam.container.Container("test1.h5", desc=name)
        cout.add("mesh",self.mesh)
        cout.add("table",self.tab)
        cout.add("data",self.data)
        cout.add("str",self.string)
        cout.write()

        cin  = steam.container.Container("test1.h5",mode="a")
        cin.read()
        cin.rm("mesh")
        cin.rm("str")
        cin.write()

        cin2 = steam.container.Container("test1.h5")
        self.assertTrue("table" in cin2.toc)
        self.assertTrue("data"  in cin2.toc)
        self.assertFalse("mesh" in cin2.toc)
        self.assertFalse("str"  in cin2.toc)
        cin2.read()
        self.assertTrue(cin2.obj["table"] == self.tab)
        self.assertTrue(cin2.obj["data" ].data.equals(self.data.data))
        self.assertFalse("mesh" in cin2.obj)
        self.assertFalse("str"  in cin2.obj)


    def test_delete(self):
        """ Make a new container, write a bunch to it, delete it all."""

        with steam.container.Container("test1.h5") as c:
            c.add("mesh",self.mesh)
            c.add("table",self.tab)
            c.add("data",self.data)
            c.add("str",self.string)
            for var in self.bad_chars:
                c.add(var,var)

        #! One thing inside of another
        with steam.container.Container("test1.h5",'a') as c:
            c.rm("mesh")
            c.rm("table")
            c.rm("data")
            c.rm("str")
            for var in self.bad_chars:
                c.rm(var)

        #! Make sure that it's empty
        ### NOTE: Since adding 'steam_version' to all containers, an "empty" 
        ###       container will have 1 item in its TOC
        with steam.container.Container("test1.h5") as c:
            self.assertEqual(len(c.toc),1)
            self.assertIn( "steam_version", c.toc )

    def test_write_exception(self):
        """ If I try and write to a read-only container, it should raise."""

        with steam.container.Container("test1.h5") as c:
            c.add("mesh",self.mesh)
            c.add("table",self.tab)

        cont = steam.container.Container("test1.h5")
        cont.read()
        cont.add("data",self.data)
        self.assertRaises(OSError,cont.write)

    def test_add_exception(self):
        """ If I try and add an object in side of another object, it should raise."""

        #! One thing inside of another
        with steam.container.Container("test1.h5") as c:
            c.add("mesh",self.mesh)
            self.assertRaises(KeyError,c.add,"mesh/table",self.tab)

        
        #! Lots of duplicates
        with steam.container.Container("test1.h5") as c:
            for var in self.bad_chars:
                c.add(var,var)
                self.assertRaises(KeyError,c.add,var,var)

    def test_string_to_hdf5( self ):
        """ Confirm we can write to HDF5 with weird column names """

        temp_cont = 'test_out.h5'
        self.files_to_remove.append( temp_cont )

        ### Create a solution and name the variable 'rho(kg/m3)'
        bad_var_name = 'rho(kg/m3)'

        ### We need another table to store additional rows
        added_tab = deepcopy( self.data )
        added_tab.data[bad_var_name] = added_tab.data['RHO']
        with steam.container.Container( temp_cont, 'w' ) as cont1:
            cont1.add( 'added_table', added_tab )
            cont1.add( 'db', self.data )

        ### Read the container that was just written
        with steam.container.Container( temp_cont, 'r' ) as cont2:
            cont2.read()

        self.assertTrue( added_tab.data.equals(cont2.obj['added_table'].data) )
        self.assertTrue( self.data.data.equals(cont2.obj['db'].data) )

    def test_selective_read(self):
        """ Read only part of a container"""

        test_cont = 'test_out.h5'
        self.files_to_remove.append( test_cont )

        with steam.container.Container( test_cont ) as c:
            c.add("mesh",self.mesh)
            c.add("table",self.tab)

        with steam.container.Container( test_cont, 'r' ) as c2:
            c2.read( ['mesh'] )

        self.assertEqual( self.mesh, c2.obj['mesh'] )

#! Run the tests
if __name__ == '__main__':
    ut.main()

