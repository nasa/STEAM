import steam

def onDiskExample():
    # CODE_START
    db = steam.database.Database("database/sol.cdat")
    db.onDisk = True
    cont = steam.container.Container("onDiskDB.h5")
    cont.add("database",db)
    cont.write()
    # Now we can load the data straight to disk
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()
    # CODE_END

def scalar_read():
    # CODE_START
    db = steam.database.Database('table/aero.cdat')
    print(db)
    print(db.data)
    # CODE_END

    OUTPUT="""
        MACH  ALPHA    RHO    QBAR    VAR1    VAR2
    0    0.0   -5.0  0.510  0.0000 -5.0000  5.0000
    1    0.0    0.0  0.510  0.0000  0.0000  0.0000
    2    0.0    5.0  0.510  0.0000  5.0000 -5.0000
    3    0.2   -5.0  0.485  0.0194 -4.9806  5.0194
    ..   ...    ...    ...     ...     ...     ...
    59   3.8    5.0  0.035  0.5054  5.5054 -4.4946
    60   4.0   -5.0  0.010  0.1600 -4.8400  5.1600
    61   4.0    0.0  0.010  0.1600  0.1600  0.1600
    62   4.0    5.0  0.010  0.1600  5.1600 -4.8400
    
    [63 rows x 6 columns]
    """ # END

    # CODE_START2
    print(db.data.columns)         # The variables
    print(db.data.MACH.min())      # Minimum MACH
    print(db.data.MACH.max())      # Maximum MACH
    np_mach = db.data.MACH.values  # MACH as NumPy array
    # CODE_END2

def scalar_set_vars():
    db = steam.database.Database('table/aero.cdat')
    # CODE_START
    db.set_indep(['MACH','ALPHA'])
    db.set_dep(['VAR1','VAR2'])
    print(db)
    # CODE_END

    OUTPUT="""
    Database Object:
     Independent Variables: ['MACH', 'ALPHA']
     Dependent Variables  : ['VAR1', 'VAR2']
     DataFrame            : (63, 6)
     Meshes               : 0 meshes
     Solution Maps        : 0 solutions
     HDF5 File            : None
     Work on Disk         : False
     Interpolation Spaces : 0 spaces
    """ # END

def scalar_mk_space():
    db = steam.database.Database('table/aero.cdat')
    db.set_indep(['MACH','ALPHA'])
    db.set_dep(['VAR1','VAR2'])
    # CODE_START
    db.new_intspace()
    # CODE_END

    OUTPUT="""
     - Scaling factors:
        - MACH: 1.0
        - ALPHA: 1.0
     - Built 2-D interpolation space
    """ # END

def scalar_interpolate():
    db = steam.database.Database('table/aero.cdat')
    db.set_indep(['MACH','ALPHA'])
    db.set_dep(['VAR1','VAR2'])
    db.new_intspace()
    # CODE_START
    case = {'MACH':3,'ALPHA':2.5,'FOO':10}
    output = db.interpolate(case)
    print(output[0].data)
    # CODE_END

    OUTPUT="""
    ALPHA     2.500
    FOO      10.000
    MACH      3.000
    VAR1      3.715
    VAR2     -1.285
    Name: 46, dtype: float64
    """ # END

    # CODE2_START
    case = {'MACH':3,'ALPHA':2.5,'FOO':10}
    output = db.interpolate(case,cols=['RHO','QBAR','VAR1','VAR2'])
    print(output[0].data)
    # CODE2_END

    OUTPUT2="""
    ALPHA     2.500
    FOO      10.000
    MACH      3.000
    QBAR      1.215
    RHO       0.135
    VAR1      3.715
    VAR2     -1.285
    Name: 46, dtype: float64
    """ # END2

    # CODE3_START
    print(output[0].get_stencil())
    # CODE3_END

    OUTPUT3="""
    i  Ind     Weight        MACH      ALPHA
    -  ---     TARGET       3.000000   2.500000
    1    46   5.00000e-01   3.000000   0.000000
    2    47   5.00000e-01   3.000000   5.000000
    """ # END3


def solution_read():
    #CODE_START
    # First, read the table of flight conditions and solution paths/types
    db = steam.database.Database('database/sol.cdat')
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()  # This actually loads the solutions into memory
    #CODE_END
    OUTPUT="""
    Read mesh: database/grid.tri
    Loading Solns : |████████████████████████████████████████████████████| 58/58 100.0%
    """ # END

def solution_mk_space():
    #CODE_START
    db = steam.database.Database('database/sol.cdat')
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()  # This actually loads the solutions into memory
    db.set_indep(['MACH','ALPHA','BETA'])
    db.new_intspace()
    traj = steam.table.Table('database/traj.cdat')
    output = db.interpolate(traj)
    #CODE_END

    OUTPUT="""
    WARNING  steam.database      : Case outside of database!  Returning no data!
    WARNING  steam.database      : Conds: 'ALPHA'=-10.000 'MACH'=2.000
    """ # END

    #CODE_START2
    print(output[0])
    #CODE_END2
    OUTPUT2="""
    Database Point:
     Data:
       MACH            : 0.1
       ALPHA           : 0.0
       BETA            : 0.0
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
    """ # END2

    #CODE_START3
    print(output[9])
    #CODE_END3
    OUTPUT3="""
    Database Point:
     Data:
     Soln:
       None
     Mesh:
       None
    """ # END3

def modify_soln_mesh():
    db = steam.database.Database('database/sol.cdat')
    db.read_mesh('database/grid.tri','TRIQ')
    db.load_soln()  # This actually loads the solutions into memory
    db.set_indep(['MACH','ALPHA','BETA'])
    db.new_intspace()
    traj = steam.table.Table('database/traj.cdat')
    output = db.interpolate(traj)
    #CODE_START
    solnid = db.data.solnid[4]  # Get the solnid for database point 11 (0-based)
    soln = db.get_soln(solnid)
    soln.data['q1'] = 2.0
    db.update_solnidx(4,soln)
    #CODE_END
    #CODE_START2
    solnid = db.data.solnid[3] 
    db.update_solnid(solnid,soln)
    #CODE_END2
    #CODE_START3
    mesh = db.get_mesh('g000000')
    mesh.transform(['t2.0x'])
    db.update_meshid('g000000',mesh)
    #CODE_END3


#scalar_mk_space()
#scalar_interpolate()
#solution_read()
#solution_mk_space()
#modify_soln_mesh()


# Call these with something like:
#import container
#func_list = [func for func in dir(container) if callable(getattr(container, func))]
#for func in func_list:
#    call = getattr(container,func)
#    call()  # Check for exception
