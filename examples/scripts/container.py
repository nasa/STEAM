import steam

# The templates in the users guide depend on the code that you
# want shown being between lines that contain:
#           CODE_START
#           CODE_END
#
# If you want output, it is probably best to contain it inside
# of a string variable named OUTPUT that is ended with an END.
# This ensures that it conforms to python syntax, does not contain
# '#' and can be cleanly parsed.
#           OUTPUT="""
#           """ END









def create(clean=True):
    # CODE_START
    # Start with a collection of objects
    tab  = steam.table.Table('table/aero.cdat')
    data = steam.database.Database('table/aero.cdat')
    mesh = steam.mesh.mk_cube(1.0,5)
    name = "Important Engineering Database"

    # Open the container (in memory) and add objects to it
    cont = steam.container.Container('test.h5',desc=name)
    cont.add('table',tab)
    cont.add('database',data)
    cont.add('meshes/mesh1',mesh)
    cont.add('meshes/mesh2',mesh)

    # Sync changes from the container to disk
    cont.write()
    print(cont)
    # CODE_END

    OUTPUT="""
    Container:
     Important Engineering Database
     HDF Path: test.h5
     HDF Mode: w
     TOC:
      database                       (toc_2 , <class 'steam.database.Database'>)
      meshes/mesh1                   (toc_3 , <class 'steam.mesh.Mesh'>)
      meshes/mesh2                   (toc_4 , <class 'steam.mesh.Mesh'>)
      table                          (toc_1 , <class 'steam.table.Table'>)
    """ # END
    if clean:
        import os
        os.remove('test.h5')

def create_context():
    # Start with a collection of objects
    tab  = steam.table.Table('table/aero.cdat')
    data = steam.database.Database('table/aero.cdat')
    mesh = steam.mesh.mk_cube(1.0,5)
    name = "Important Engineering Database"

    # CODE_START
    # Open the container (in memory) and add objects to it
    with steam.container.Container('test.h5',desc=name) as cont:
      cont.add('name',name)
      cont.add('table',tab)
      cont.add('database',data)
      cont.add('meshes/mesh1',mesh)
      cont.add('meshes/mesh2',mesh)
    # CODE_END

    # OUTPUT
    # END
    import os
    os.remove('test.h5')

def read_container():
    create(False)
    # CODE_START
    # Open container made previously
    with steam.container.Container('test.h5') as cont:
        cont.read() # Load the contents into memory
        table = cont.obj['table']
        mesh  = cont.obj['meshes/mesh1']
    print(mesh)
    # CODE_END

    OUTPUT="""
    Reading HDF5 File 'test.h5'
    Python Mesh Object:
       Points    : 98
       Elements  : 192
       # of Comps: 6
       Int Components: 1 2 3 4 5 6
       Str Components:
    """ #END
    import os
    os.remove('test.h5')

def view_contents():
    create(False)
    # CODE_START
    # Open container made previously
    cont = steam.container.Container('test.h5')
    print(cont)
    # CODE_END

    OUTPUT="""
    Container:
     Important Engineering Database
     HDF Path: test.h5
     HDF Mode: r
     TOC:
      database                       (toc_2 , <class 'steam.database.Database'>)
      meshes/mesh1                   (toc_3 , <class 'steam.mesh.Mesh'>)
      meshes/mesh2                   (toc_4 , <class 'steam.mesh.Mesh'>)
      table                          (toc_1 , <class 'steam.table.Table'>)
    """ #END
    import os
    os.remove('test.h5')

def remove_items():
    create(False)
    # CODE_START
    # Delete items from previous container
    with steam.container.Container('test.h5',mode='a') as cont:
        cont.read() # Load the contents into memory
        cont.rm('database')
        cont.rm('meshes/mesh1')
    # Reopen and confirm reduced contents
    with steam.container.Container('test.h5') as cont:
        print(cont)
    # CODE_END

    OUTPUT="""
    Reading HDF5 File 'test.h5'
    Container:
     Important Engineering Database
     HDF Path: test.h5
     HDF Mode: r
     TOC:
      table                          (toc_1 , <class 'steam.table.Table'>)
      meshes/mesh2                   (toc_4 , <class 'steam.mesh.Mesh'>)
    """ #END
    import os
    os.remove('test.h5')


#create()
#create_context()
#read_container()
#view_contents()
#remove_items()
