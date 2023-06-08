#!/usr/bin/env python
""" Create random dataframes for use in test_table.py.

This function can be used at any time to regenerate the files necessary for
testing the `steam.table.Table` class.
"""

import pandas as pd
import random

from code import interact

### Define range of each variable
var_limits = {
    "Time(sec)"   : (1, 1000),
    "Alt(km)"     : (20, 60),
    "Vel(m/s)"    : (3e3, 7.5e3),
    "Press(Pa)"   : (20e3, 40e3),
    "Temp(K)"     : (200, 300),
    "Rho(Kg/m3)"  : (5e-4,1e-3),
    "Alpha(Deg)"  : (10, 40),
    "Beta(Deg)"   : (-10, 10),
    "AlphaT(Deg)" : (10, 50),
    "Phi(Deg)"    : (0, 180),
}
    
#######################################

def mk_dataframe( traj_names ) -> pd.DataFrame:
    """ Given a list containing names of the "trajectories" in the table of
    fake data that we're creating, build a pandas.DataFrame with a point
    for each trajectory.
    """
    
    table_dict = {'Traj' : traj_names}

    n_traj = len(traj_names)
    for var, limits in var_limits.items():
        ### Number of entries for each variable must equal number of trajs
        table_dict[ var ] = [ random.uniform(*limits) for i in range(n_traj) ]

    return pd.DataFrame.from_dict( table_dict )

#######################################
def main():

    ### Define output file names and their fake trajectory names
    file_dict = { 'table1.dat' : ['Tr01','Tr02','Tr03','Tr04','Tr05','Tr06'],
                  'table2.dat' : ['Tr07','Tr08','Tr09'],
                  'table2_bad.dat' : ['Tr10','Tr11','Tr12'],
                }

    for f_name, t_names in file_dict.items():
        df = mk_dataframe( t_names )

        if 'bad' in f_name:
            ### Deliberately creating a file that can't be added to table1.dat
#            interact( local = dict( globals(), **locals() ) )
            df.drop( columns = 'Vel(m/s)', inplace=True )

        df.to_csv( f_name, sep='\t', float_format='%13.5f', index=False )
#######################################

if __name__ == "__main__":
    main()
