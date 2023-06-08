import steam

def extract_point():
    # CODE_START
    db = steam.database.Database('database/sol.cdat')
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()  # This actually loads the solutions into memory
    ind0 = db.get_point(0)
    print(ind0)
    # CODE_END

    OUTPUT="""
    Database Point:
     Data:
       MACH            : 0.0
       ALPHA           : -5.0
       BETA            : -5.0
       linear          : -2.5
       quad            : 12.5
     Soln:
          Data Size       : 6
          Variables       : q1 q2 q3 q4
          Node-based
          Mesh             :
             Points    : 6
             Elements  : 4
             # of Comps: 1
             Int Components: 1
             Str Components:
             No Static Solution.
     Mesh:
          Points    : 6
          Elements  : 4
          # of Comps: 1
          Int Components: 1
          Str Components:
          No Static Solution.
    """ # END

def expand_alpha():
    db = steam.database.Database('database/sol.cdat')
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()  # This actually loads the solutions into memory
    ind0 = db.get_point(0)
    # CODE_START
    ind0.data.ALPHA = -180
    db.add_point(ind0)
    print(db.data.iloc[-1])
    # CODE_END

    OUTPUT="""
    MACH               0
    ALPHA           -180
    BETA              -5
    linear          -2.5
    quad            12.5
    soln_type        NaN
    soln_path        NaN
    meshid       g000000
    solnid       s000000
    Name: 58, dtype: object
    """ # END

    db = steam.database.Database('database/sol.cdat')
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()  # This actually loads the solutions into memory
    ind0 = db.get_point(0)
    # CODE_START2
    cases = db.data[db.data["ALPHA"] == -5]
    for (i,case) in cases.iterrows():
      new_pt = db.get_point(i)
      new_pt.data.ALPHA = -180.0
      db.add_point(new_pt)
    # CODE_END2
    print(db)

    OUTPUT2 = """
    Database Object:
     Independent Variables: []
     Dependent Variables  : []
     DataFrame            : (62, 9)
     Meshes               : 1 meshes
        g000000 :       6 points
     Solution Maps        : 58 solutions
        s000023 : 6 points,   4 vars
        ...
        s000057 : 6 points,   4 vars
     HDF5 File            : None
     Work on Disk         : False
     Interpolation Spaces : 0 spaces
     """ # END2

def update_mesh():
    db = steam.database.Database('database/sol.cdat')
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()  # This actually loads the solutions into memory
    # CODE_START
    print(db.data.meshid)
    for (i,point) in enumerate(db):
        point.mesh.transform(["t{}x".format(i)])
        point.update_mesh()
    print(db.data.meshid)
    # CODE_END

    OUTPUT = """
    0     g000000
    1     g000000
    ...     ...
    56    g000000
    57    g000000
    
    Translating in x by 0.0 units
    Translating in x by 1.0 units
                ...
    Translating in x by 56.0 units
    Translating in x by 57.0 units
    
    0     g000000
    1     g000001
    ...     ...
    56    g000056
    57    g000057
     """ # END

#extract_point()
#expand_alpha()
#update_mesh()


# Call these with something like:
#import container
#func_list = [func for func in dir(container) if callable(getattr(container, func))]
#for func in func_list:
#    call = getattr(container,func)
#    call()  # Check for exception

#extract_point()
#expand_alpha()


# Call these with something like:
#import container
#func_list = [func for func in dir(container) if callable(getattr(container, func))]
#for func in func_list:
#    call = getattr(container,func)
#    call()  # Check for exception
