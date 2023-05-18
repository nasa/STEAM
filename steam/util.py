# -*- coding: utf-8 -*-
""" Utility Module

Module of common utilities for STEAM.
"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import numpy  as np
from os import popen

def hasher(in_string):
    """ Use a stable hasher to return the hex hash of
    the input string.  The returned string should be
    unique enough for stable comparisons.

    xxhash is the faster option (like 25x faster).

    Args:
        in_string (:obj:`bytes`): Byte string to hash.

    Returns:
        :obj:`str`: Hash of the input string.
    """

    out = None;

    #import xxhash
    #h=xxhash.xxh64()
    ##print(type(in_string))
    ##print(in_string.flags)
    #h.update(in_string)
    #out = h.hexdigest()
    #print(out)
    #import hashlib
    #print(hashlib.sha1(in_string).hexdigest())

    try:
        import xxhash
        h=xxhash.xxh64()
        h.update(in_string)
        out = h.hexdigest()
    except:
        raise Exception(" ** Could not form hash using xxhash!")
    #try:
    #    import hashlib
    #    out = hashlib.sha1(in_string).hexdigest()
    #except:
    #    raise Exception(" ** Could not form hash using hashlib!")

    return out


def vars_in_list(vars,cols):
    """Check to make sure that the variables exists in
    the provided list

    Args:
        cols (:obj:`list`): List of things to look in.
        vars (:obj:`list`): List of variables to ensure exist.

    Returns:
        :obj:`bool`: All vars in cols? (T/F)
    """

    stat = True
    for i in vars:
        if (not i in cols):
            stat = False

    return stat

def get_var_cap(varlist,var):
    """Look in the items in the varlist to see if var exists.  If
    it doesn't, then assume that it does and that capitalization
    is the issue.  If that fails, then raise an error.  Success
    returns the proper capitalization of the variable.

    Args:
     varlist (:obj:`list`): List of variable names to check.
     var     (:obj:`str`): Variable to find in list.

    Returns:
     varcap (:obj:`str`): Properly capitalized variable string.
    """

    varcap = var

    if not var in varlist:
        ### Now check against different capitalization
        lowlist = list(map(str.lower, varlist))
        lowvar  = var.lower()

        try:
            varind = lowlist.index(lowvar)
        except:
            raise IOError(
                "Variable {} not in varlist!\n".format(var)
                )
        varcap = varlist[varind]

    return varcap

def run_command(cmd):
    """Runs the given command locally and returns the output, err and exit_code. 

      From: http://stackoverflow.com/questions/7389662/link-several-popen-commands-with-pipes

    Args:
        cmd (:obj:`str`): Command to run
    Returns:
        (:obj:`str`): Standard Output
        (:obj:`str`): Standard Error
        (:obj:`int`): Exit Code
    """

    from subprocess import Popen, PIPE
    import shlex
    if "|" in cmd:
       cmd_parts = cmd.split('|')
    else:
       cmd_parts = []
       cmd_parts.append(cmd)
    i = 0
    p = {}
    for cmd_part in cmd_parts:
       cmd_part = cmd_part.strip()
       if i == 0:
          p[i]=Popen(shlex.split(cmd_part),stdin=None, stdout=PIPE, stderr=PIPE)
       else:
          p[i]=Popen(shlex.split(cmd_part),stdin=p[i-1].stdout, stdout=PIPE, stderr=PIPE)
       i = i +1
    (output, err) = p[i-1].communicate()
    exit_code = p[0].wait()
    return str(output), str(err), exit_code



def progress_bar (iteration, total, prefix = '', suffix = '', decimals = 1, length = None, fill = '█', frac = 10):
    """ Call in a loop to create terminal progress bar.

    The "iteration" input should start with 1 rather than 0; otherwise the bar
    will never reach 100%.

    Adapted from here: http://stackoverflow.com/a/34325723

    Example Usage::

        # Initial call to print 0% progress
        i= 0; l = len(items);
        progress_bar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)
        for item in items:
            # Do stuff...
            # Update Progress Bar
            i += 1
            progress_bar(i, l, prefix = 'Progress:', suffix = 'Complete', length = 50)

    Args:
        iteration   (:obj:`int`)  : current iteration
        total       (:obj:`int`): total iterations
        prefix      (:obj:`str`, Optional): prefix string
        suffix      (:obj:`str`, Optional): suffix string
        decimals    (:obj:`int`, Optional): positive number of decimals in percent complete
        length      (:obj:`int`, Optional): character length of bar
        fill        (:obj:`str`, Optional): bar fill character
        frac        (:obj:`int`, Optional): Fraction of the bar to print (number of steps)
                                            -1 is every iteration
    """

    # We're only going to write things in 10% increments.  That is becuase, it turns out, this can
    # be expensive!
    pct = max(int(total / frac),1)
    if (iteration % pct != 0 and iteration != total and iteration != 1):
        return

    if (iteration > total):
        return

    ### Get window width for sizing the bar
    try:
        #_, width = popen( 'stty size', 'r').read().split()
        with popen( 'stty size', 'r' ) as f:
            _, width = f.read().split()
    except:
        width = 84
    width = int(width)

    ### Determine how long the bar should be to keep output on one line
    if length is None:
        length = width - (12 + len(prefix) + len(suffix) + 
                        2*len(str(total)) + decimals )
    length = 0 if length < 0 else length    ### length must be >= 0

    n_percent = 4 + len(str(decimals))
    percent = ( "{0:" + str(n_percent) + '.' + str(decimals) + "f}").format(
                100 * (iteration / float(total)) )
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    out_string = '\r%s |%s| %' + str(len(str(total))) + 'd/%s %s%% %s'
    print( out_string % (prefix, bar, iteration, total, percent, suffix), 
           end = '\r')

    # Print New Line on Complete
    if iteration == total: 
        print()



#! ## HDF5 functions

def write_hdf5_pnode(store,path,obj,compress=False):
    """Write a pickled node to an HDF5 file.
    
    Args:
        store (:obj:`HDFStore`): HDF5 Store to put it in
        path  (:obj:`str`):   Path to where h5 file to put it.
        obj   (:obj:`object`):   Object to pickle and store
        compress  (:obj:`Boolean`):   Compress? [F]
    
    """

    import pickle
    import tables
    from tables.nodes import filenode

    #! I need a path since I can't just use the pandas HDF5store
    file = store.filename

    import os.path
    import pathlib

    path = path.replace("//","/")
    dir_name = os.path.dirname( path)
    var_name = os.path.basename(path)

    #! Open a filenode at /toc
    h5file = tables.open_file(file,'a')

    #! Long story short, I need to make group that don't exist from the root path
    #! all the way to the object that I want to write.
    p = pathlib.PurePath(dir_name)
    comb_path = p.parts[0]
    for part in p.parts[1:]:
        exists = store.get_node(comb_path+"/"+part)
        #print(part,'->',comb_path+"/"+part,'->',exists)
        if not exists:
            h5file.create_group(comb_path,part)
        #print(part,'->',comb_path+"/"+part,'->',exists)
        comb_path = os.path.join(comb_path,part)

    try:
        fnode = filenode.new_node(h5file, where=dir_name, name=var_name)
    except tables.exceptions.NodeError:
        h5file.remove_node(path)
        fnode = filenode.new_node(h5file, where=dir_name, name=var_name)

    if compress:
        import zlib
        fnode.write(zlib.compress(pickle.dumps(obj)))
    else:
        #! Pickle and dump it
        pickle.dump(obj,fnode)

    #! Close files
    fnode.close()
    h5file.close()


def read_hdf5_pnode(store,path,compress=False):
    """Read a pickled node from an HDF5 file.
    
    Args:
        store (:obj:`HDFStore`): HDF5 Store to read from
        path  (:obj:`str`):   Path in h5 file
        compress  (:obj:`Boolean`):   Compressed?
    
    Returns:
        (:obj:`object`):   Unpickled object
    """

    import pickle
    import tables
    from tables.nodes import filenode

    import os.path

    path = path.replace("//","/")

    #! Open a filenode at /toc
    node   = store.get_node(path)
    fnode = filenode.open_node(node, 'r')

    if compress:
        import zlib
        obj = pickle.loads(zlib.decompress(fnode.read()))
    else:
        #! Read the object
        obj = pickle.load(fnode)

    #! Close files
    fnode.close()

    return obj

def hdf5_key_to_string(key):
    """ Convert HDF5 key to string. 
    
    Args:
        key (:obj:`str`): HDF5 key
    
    Returns:
        (:obj:`str`): Original string
    """

    import re
    string = key
    #! Has to begin with a letter or _
    string = re.sub('^_h5s_','',string)
    #! Cannot have periods
    string = re.sub('_p5s_','.',string)
    #! Cannot have dashes
    string = re.sub('_d5s_','-',key)
    #! Cannot have plus
    string = re.sub('_pl5s_','+',key)

    return string

def string_to_hdf5_key(string):
    """ Convert string to valid HDF5 key. 
    
    This is necessary since HDF5 keys have to begin with a letter
    and cannot have periods or '-'.  So we do some substitution here.

    Args:
        string (:obj:`str`): Original string
    
    Returns:
        (:obj:`str`): HDF5 key
    """

    import re
    key    = string
    #! Has to begin with a letter or _
    #! Can allow a '/' in case it's a path
    key    = re.sub('^([^a-zA-Z_/])',r'_h5s_\1',key)
    #! Cannot have periods
    key    = re.sub('\.','_p5s_',key)
    #! Cannot have dashes
    key    = re.sub('\-','_d5s_',key)
    #! Cannot have dashes
    key    = re.sub('\+','_pl5s_',key)

    return key


def df_to_cdat(filename,df,header=None,var_pre="#",ovars=None):
    """ Write a Pandas Dataframe (neatly) to a cdat file.
    
    By default this makes a columnized data file with the variable names in the
    first line prefaced by a '#'.  To add additional comment lines above the
    variable names, use the header keyword.

    To write a tecplot ASCII file, pass var_pre='variables = ' and provided the
    variable/column names do not contain spaces it should be readable as a .dat file.

    Args:
        filename (:obj:`string`): File to write.
        df (:obj:`~pandas.DataFrame`): dataframe to write
        header (:obj:`string`) : String to write before the variable names and data. [""]
        var_pre (:obj:`string`) : Prefix to place on the line containing variable names. [#]
        ovars (:obj:`list`): Optional, list/order of variables to write.
        
    """

    if ovars is None:
        ovars = list(df.columns)

    #! Minimum length of variable (based on name)
    min_len = dict()
    for (key,hlen) in zip(ovars,[len(name) for name in ovars]):
        min_len[key] = hlen

    #! This is a bit, just a bit of a hack but I want to identify which
    #! columbs are floats, ints, or strings ('objects' in pandas)
    hform = dict() # Header format
    dform = dict() # Data   format

    #! Floats
    for var in df.columns[df.dtypes == float]:
        if var not in ovars:
            continue
#       print("Float Var = {}".format(var))
        hlen = max(16,min_len[var])
        dform[var] = '%{}.12g'.format(hlen)
        hform[var] = '^{}s'.format(hlen)

    #! Ints
    for var in df.columns[df.dtypes == int]:
        if var not in ovars:
            continue
#       print("Int Var = {}".format(var))
        hlen = max(16,min_len[var])
        dform[var] = '%{}d'.format(hlen)
        hform[var] = '^{}s'.format(hlen)

    #! Strings
    for var in df.columns[df.dtypes == object]:
        if var not in ovars:
            continue
#       print("String Var = {}".format(var))
        lens = [len(str(item)) for item in df[var]]
        hlen = max(max(lens),min_len[var])
        dform[var] = '%{}s'.format(hlen)
        hform[var] = '^{}s'.format(hlen)

    ohead = []
    oform = []
    for var in ovars:
        ohead.append(("{:"+"{}".format(hform[var])+"}").format(var))
        oform.append(dform[var])

    ohead = " ".join(ohead)
    ohead = var_pre+" "+ohead
    
    if header is not None:
        ohead = header+"\n"+ohead

    #! The call sets 'comments=''' so that numpy does not prepend everything with '#'
    np.savetxt(filename,df[ovars].values,fmt=oform,header=ohead,comments='')

# Logging Utilities

def timer(logger, level=None):
    """ Simple decorator to return the timing for a function.
        
        Usage should just be to add @steam.util.timer(logger) before a
        function declaration
    """
    import time
    from functools import wraps

    if level is None:
        level = logging.DEBUG

    def wrap(fn):
        @wraps(fn)
        def decorator(*args, **kwargs):
            start   = time.time()
            logger.log(level, "Starting '{}'".format(
                            fn.__name__,
                            )
                       )
            result  = fn(*args, **kwargs)
            duration = time.time() - start
            logger.log(level, "'{}' took {:12f} ms".format(
                            fn.__name__,
                            duration * 1000
                            )
                       )
            return result
        return decorator

    return wrap

def read_docs():
    """ Open the STEAM documentation in a web browser.

    This function uses the webbrowser module to launch a browser tab with the
    STEAM documentation.  There are several reasons that this function will
    fail, including if the documentation is not compiled, or if the user's
    Python distribution does not have the webbrowser module.
    """

    import steam
    import webbrowser as wb

    ### Find the location of the documentation, based on the location of the
    ### STEAM module
    doc_index = steam.__path__[0] + '/../docs/build/html/index.html'

    return wb.open( doc_index, new=2 )
