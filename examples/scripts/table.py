import steam

def load():
    # CODE_START
    table = steam.table.Table('table/aero.cdat')
    print(table.data[0:3])
    # CODE_END

    OUTPUT="""
       MACH  ALPHA   RHO  QBAR  VAR1  VAR2
    0   0.0   -5.0  0.51   0.0  -5.0   5.0
    1   0.0    0.0  0.51   0.0   0.0   0.0
    2   0.0    5.0  0.51   0.0   5.0  -5.0
    """ # END

def min_max():
    table = steam.table.Table('table/aero.cdat')

    # CODE_START
    table.data.min()
    table.data.max()
    # CODE_END

    OUTPUT="""
    ALPHA   -5.00
    RHO      0.01
    QBAR     0.00
    VAR1    -5.00
    VAR2    -5.00
    dtype: float64

    MACH     4.0000
    ALPHA    5.0000
    RHO      0.5100
    QBAR     1.2544
    VAR1     6.2544
    VAR2     6.2544
    dtype: float64
    """ # END
    
def rm_column():
    table = steam.table.Table('table/aero.cdat')

    # CODE_START
    table.remove_vars(['QBAR','RHO'])
    print(list(table.data.columns))
    # CODE_END

    OUTPUT="""
    ['MACH', 'ALPHA', 'VAR1', 'VAR2']
    """ # END
   
def add_column():
    table = steam.table.Table('table/aero.cdat')
    table.remove_vars(['QBAR','RHO'])

    # CODE_START
    import numpy as np
    new      = np.ones((table.data.shape[0],2))
    new[:,1] = 2.0
    table.add_vars(["ONE","TWO"],new)
    print(table.data[0:3])
    # CODE_END

    OUTPUT="""
       ALPHA  MACH  ONE  TWO  VAR1  VAR2
    0   -5.0   0.0  1.0  2.0  -5.0   5.0
    1    0.0   0.0  1.0  2.0   0.0   0.0
    2    5.0   0.0  1.0  2.0   5.0  -5.0
    """ # END
 
#rm_column()
#add_column()


# Call these with something like:
#import container
#func_list = [func for func in dir(container) if callable(getattr(container, func))]
#for func in func_list:
#    call = getattr(container,func)
#    call()  # Check for exception
