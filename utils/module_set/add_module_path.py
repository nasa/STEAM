#!/usr/bin/env python
# add_module_path.py
""" Provide the user with a line for their .bashrc or .cshrc file to enable
use of the provided STEAM module.

Note that the provided module will only work on systems that use the Lmod
module system.  This script, being a relatively simple Python routine, will
work anywhere, but it should only be executed in an Lmod environment.

    Args:
        rcpath (Optional): Path to the user's bashrc/cshrc/etc file.  The "module use" command will be appended to this file in order to enable use of the provided STEAM module in every shell.
"""

import os

##################
def get_command():
    """ Return a string containing the appropriate "module use" command
    """
    ### We want the absolute path to the directory in which add_module_path.py
    ###     lives, regardless of the directory from which it's called.
    module_dir = os.path.dirname( os.path.abspath(__file__) )
    mod_command = 'module use {}'.format( module_dir )
    return mod_command
##################

if __name__ == '__main__':
    ### Establish the argparser to read command line args
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( "--rcpath", help="""Path to the user's """
             + """bashrc/cshrc/etc file.  The "module use" command will be """
             + "appended to this file in order to enable use of the provided "
             + "STEAM module in every shell." )
    args = parser.parse_args()

    mod_command = get_command()
    print( mod_command )

    if args.rcpath:
        ### Append mod_command to provided bashrc/cshrc path
        with open(args.rcpath, 'a') as rc_file:
            rc_file.write( '### Include path to STEAM module in module path\n')
            rc_file.write( mod_command + "\n" )
