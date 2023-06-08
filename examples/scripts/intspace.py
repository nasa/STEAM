import steam
import numpy as np
import pandas as pd

def database():
    # CODE_START
    import numpy as np
    import pandas as pd
    sphere = steam.mesh.mk_sphere (radius=0.5, npts=11)
    solns = []
    for i in range (4):
        data = np.ones (sphere.xyz_pt.shape [0]) * i
        soln = steam.solution.Solution (mesh=sphere)
        soln.add_var (var='pressure', data=data)
        solns.append (soln)
    points = [steam.database.DBPoint ( data=pd.Series ( {'Mach' : 1.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [0]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [1]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 2.0, 'alpha' : 0.0 } ), soln=solns [2]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 1.0 } ), soln=solns [3])]
    db = steam.database.Database (description='A very simple 3D database that contains a single simplex')
    for point in points:
        db.add_point (point)
    print (db)
    # CODE_END

    OUTPUT="""
    Database Object:
     A very simple 3D database that contains a single simplex
     Independent Variables: []
     Dependent Variables  : []
     DataFrame            : (4, 5)
     Meshes               : 1 meshes
        g000000 :     200 points
     Solution Maps        : 4 solutions
        s000002 : 200 points,   1 vars
        ...
        s000000 : 200 points,   1 vars
     HDF5 File            : None
     Work on Disk         : False
     Interpolation Spaces : 0 spaces
    """ # END

def new_intspace():

    # CODE_START
    import numpy as np
    import pandas as pd
    sphere = steam.mesh.mk_sphere (radius=0.5, npts=11)
    solns = []
    for i in range (4):
        data = np.ones (sphere.xyz_pt.shape [0]) * i
        soln = steam.solution.Solution (mesh=sphere)
        soln.add_var (var='pressure', data=data)
        solns.append (soln)
    points = [steam.database.DBPoint ( data=pd.Series ( {'Mach' : 1.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [0]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [1]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 2.0, 'alpha' : 0.0 } ), soln=solns [2]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 1.0 } ), soln=solns [3])]
    db = steam.database.Database (description='A very simple 3D database that contains a single simplex')
    for point in points:
        db.add_point (point)
    db.new_intspace (key='mqa', indep=['Mach', 'qbar', 'alpha'], scale='auto', opts={'qbar':"log10"})
    # CODE_END

    OUTPUT="""
     - Variable Operations:
        - qbar -> log_10(qbar)
     - Scaling factors:
        - Mach: 1.0
        - qbar: 3.321928094887362
        - alpha: 1.0
    """ # END
    
def interpolate():

    # CODE_START
    import numpy as np
    import pandas as pd
    sphere = steam.mesh.mk_sphere (radius=0.5, npts=11)
    solns = []
    for i in range (4):
        data = np.ones (sphere.xyz_pt.shape [0]) * i
        soln = steam.solution.Solution (mesh=sphere)
        soln.add_var (var='pressure', data=data)
        solns.append (soln)
    points = [steam.database.DBPoint ( data=pd.Series ( {'Mach' : 1.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [0]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [1]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 2.0, 'alpha' : 0.0 } ), soln=solns [2]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 1.0 } ), soln=solns [3])]
    db = steam.database.Database (description='A very simple 3D database that contains a single simplex')
    for point in points:
        db.add_point (point)
    db.new_intspace (key='mqa', indep=['Mach', 'qbar', 'alpha'], scale='auto', opts={'qbar':"log10"})
    case = {'Mach':1.9, 'qbar':1.1, 'alpha':0.1}
    intpoint = db.interpolate (case, space='mqa')
    print (intpoint)
    # CODE_END

    OUTPUT="""
     - Variable Operations:
        - qbar -> log_10(qbar)
     - Scaling factors:
        - Mach: 1.0
        - qbar: 3.321928094887362
        - alpha: 1.0
    [
    Database Point:
     Data:
       Mach            : 1.9
       alpha           : 0.1
       qbar            : 1.1
     Stencil:
       i  Ind     Weight        Mach       qbar      alpha   
       -  ---     TARGET       1.900000   1.100000   0.100000
       0     3   1.00000e-01   2.000000   1.000000   1.000000
       1     2   1.37504e-01   2.000000   2.000000   0.000000
       2     1   6.62496e-01   2.000000   1.000000   0.000000
       3     0   1.00000e-01   1.000000   1.000000   0.000000
     Soln:
          Data Size       : 200
          Variables       : pressure
          Node-based
          Mesh             : 
             Points    : 200
             Elements  : 378
             # of Comps: 1
             Int Components: 1
             Str Components: 
             No Static Solution.
     Mesh:
          Points    : 200
          Elements  : 378
          # of Comps: 1
          Int Components: 1
          Str Components: 
          No Static Solution.
       ]
    """ # END
   
def get_weights():

    # CODE_START
    import numpy as np
    import pandas as pd
    sphere = steam.mesh.mk_sphere (radius=0.5, npts=11)
    solns = []
    for i in range (4):
        data = np.ones (sphere.xyz_pt.shape [0]) * i
        soln = steam.solution.Solution (mesh=sphere)
        soln.add_var (var='pressure', data=data)
        solns.append (soln)
    points = [steam.database.DBPoint ( data=pd.Series ( {'Mach' : 1.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [0]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 0.0 } ), soln=solns [1]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 2.0, 'alpha' : 0.0 } ), soln=solns [2]), \
              steam.database.DBPoint ( data=pd.Series ( {'Mach' : 2.0, 'qbar' : 1.0, 'alpha' : 1.0 } ), soln=solns [3])]
    db = steam.database.Database (description='A very simple 3D database that contains a single simplex')
    for point in points:
        db.add_point (point)
    db.new_intspace (key='mqa', indep=['Mach', 'qbar', 'alpha'], scale='auto', opts={'qbar':"log10"})

    cases = {'Mach':[1.80, 1.85, 1.90], 'qbar':[1.1, 1.15, 1.20], 'alpha':[0.0, 0.05, 0.10]}
    df    = pd.DataFrame (data=cases)
    db.intspace ['mqa'].get_weights (df)
    # CODE_END

    OUTPUT="""
     - Variable Operations:
        - qbar -> log_10(qbar)
     - Scaling factors:
        - Mach: 1.0
        - qbar: 3.321928094887362
        - alpha: 1.0
    ([array([3, 2, 1, 0]), array([3, 2, 1, 0]), array([3, 2, 1, 0])], [array([ 0.        ,  0.13750352,  0.66249648,  0.2       ]), array([ 0.05      ,  0.20163386,  0.59836614,  0.15      ]), array([ 0.1       ,  0.26303441,  0.53696559,  0.1       ])])
    """ # END

#def interpolate_numpy_ndarray():
#
#    # CODE_START
#    import math
#    n = 11
#    xyz = np.zeros ((n*n*n, 3))
#    data = np.zeros ((n*n*n, 1))
#    c = 0
#    for i in range (n):
#        for j in range (n):
#            for k in range (n):
#                x = -1.0 * math.pi / 2 + (math.pi / (n - 1) * i)
#                y = -1.0 * math.pi / 2 + (math.pi / (n - 1) * j)
#                z = -1.0 * math.pi / 2 + (math.pi / (n - 1) * k)
#                dist = math.sqrt (x**2 + y**2 + z**2)
#                xyz [c][0] = x
#                xyz [c][1] = y
#                xyz [c][2] = z
#                data [c][0] = math.cos (dist)
#                c += 1
#                
#    # CODE_END
