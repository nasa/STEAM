""" Module for the Models Class.

Basically a collection of routines that can be used to manipulate
solutions.  This may not be that useful down the road but is
currently something to abstract the pandas manipulation necessary
to interface between meshes, components, solutions, and variables.
"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import numpy  as np
import steam.util

###############################################################

def scalar_to_component(soln,comp,var_to,scalar,method='repl'):
    """ Modifies data on a component in a soln dataframe by a scalar."

    Args:
        soln (:obj:`~steam.solution.Solution`): Solution to augment
        comp (:obj:`int`): Integer component number, string component name, or a list containing either or both of those, indicating which component(s) to augment
        var_to (:obj:`str`): Variable to augment
        scalar (:obj:`float`): Scalar value to use in augmentation
        method (:obj:`str`): Method for augmentation.  Options: "add", "sub", "mult", "div", and "repl".  Defaults to "repl".
    """

    #! TODO: Check to make sure that soln is element-based (??)
    if (soln.store_at.upper().find('ELEM') != 0):
        raise Exception(
              "Solution must be element based to manipulate components!" )

    #! TODO: Check to make sure that var_to exist in soln/data
    if (not var_to in soln.data.columns):
        raise ValueError("Variable {} not in solution!".format(var_to))

    ### Allow for a list of components to be given
    if isinstance(comp, str) or isinstance( comp, int ):
        comp = [comp]

    try:
        loc = soln.mesh.get_comp(comp)
    except:
        raise "Could not find component, '"+str(comp)+"', in mesh!"

    combo = 0

    #! Check for the operation and update combo
    if   (method.lower() == "add"  or method == "+"):
        combo = soln.data[var_to].loc[loc] + scalar
    elif (method.lower() == "sub"  or method == "-"):
        combo = soln.data[var_to].loc[loc] - scalar
    elif (method.lower() == "mult" or method == "*"):
        combo = soln.data[var_to].loc[loc] * scalar
    elif (method.lower() == "div"  or method == "/"):
        combo = soln.data[var_to].loc[loc] / scalar
    elif (method.lower() == "repl"):
        #! Get the shape right, then overwrite
        combo = soln.data[var_to].loc[loc]
        combo[:] = scalar
    else:
        print("* Error: method not recognized:",method)
        raise IOError

    #! Apply the update to the soln data
    soln.data.update(combo)
 
    return

#sub = foo2.iloc[[1,3,5,7]]
#combo = foo.iloc[sub.index] + sub
#foo.update(combo['A'])
#
#foo[:] = 1

###############################################################

def soln_by_soln( in_soln, soln_cols, model_soln, model_col,
                  operation = 'MULT', in_place = False, model_index=None ):
    """ Augment a solution with factors contained in another solution.

    An implicit assumption in this function is that `in_soln` and `model_soln`
    are on similar meshes.  Mesh equality is not enforced because of potential
    differences in components or other properties that would be irrelevant,
    but unexpected behavior may occur if 
    in_soln.data.shape[0] != model_soln.data.shape[0]

    The `operation` input defines how the values in the model solution will
    be used.  The `model_index` input allows the user to specify a subset
    of the solutions to augment.  This is useful when one is interested in
    modifying only part of a solution.
        "MULT"      The data in `model_soln` will be multiplied by the data in `in_soln`
        "ADD"       The data in `model_soln` will be added to the data in `in_soln`
        "REPLACE"   The data in `model_soln` will overwrite the data in `in_soln`
        "MAX"       The maximum of the data in `model_soln` and the data in `in_soln` will be used
        "MIN"       The minimum of the data in `model_soln` and the data in `in_soln` will be used

    Args:
        in_soln (:obj:`~steam.solution.Solution`): The solution to be augmented
        soln_cols (:obj:`list`): List containing the solution columns to be augmented
        model_soln (:obj:`~steam.solution.Solution`): The solution containing augmentation factors
        model_col (:obj:`str` or :obj:`int`): Column in model_soln containing model factors.  Can be a string containing the column name or an integer column index.
        operation (:obj:`str`): Flag indicating how augmentation factors should be applied.  Defaults to 'MULT' for multiplicative augmentation.
        in_place (:obj:`bool`): Flag indicating whether to modify the input solution in place. Defaults to False.
        model_index (:obj:`list`): Index of items in the solutions subject to augmentation.  Defaults to all items.
    Returns:
        (:obj:`~steam.solution.Solution`): A solution containing the augmented/manipulated data
    """

    ### Initialize the output solution
    if in_place:
        out_soln = in_soln
    else:
        out_soln = steam.solution.Solution( copy = in_soln )

    if model_index:
        ### An index has been provided; no action required
        pass
    else:
        ### The model_soln and in_soln should have identical indices
        model_index = model_soln.data.index

    ### Get model solution column name, not index
    if isinstance( model_col, int ):
        model_col = model_soln.data.columns[ model_col ]

    ### Get the data used to augment the in_soln at the appropriate indices
    af_data = model_soln.data[model_col].loc[model_index]

    logger.debug('Augmenting solution by column:  {}'.format(af_data.name))

    ### Define the method of augmentation
    operation = operation.upper()
    if operation == 'MULT':
        operator = lambda x,y: x * y
    elif operation == 'ADD':
        operator = lambda x,y: x + y
    elif operation == 'REPLACE':
        operator = lambda x,y: y
    elif operation == 'MAX':
        operator = lambda x,y: max(x,y)
    elif operation == 'MIN':
        operator = lambda x,y: min(x,y)
    else:
        raise ValueError('operation input must be "MULT", "ADD", or "REPLACE"')

    ### Augment the appropriate columns in the solution
    for col in soln_cols:
        out_soln.data[col].loc[model_index] = operator( 
                                out_soln.data[col].loc[model_index], af_data )
        logger.debug('Augmenting solution column:  {}'.format( col ) )

    return out_soln

###############################################################
