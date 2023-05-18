""" Module for the Interpolation
"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import pandas as pd
import numpy  as np
import steam
from   scipy.interpolate import LinearNDInterpolator
from   scipy.spatial     import Delaunay
from   scipy.spatial     import cKDTree

#from   code              import interact
#interact( local = dict( globals(), **locals() ) )

class Transform():
    """Class definition for Transformations.

    Transformations are indicies/weights for converting between
    one grid and another.
    """

    def __init__(self,source=None,target=None,
                 elements=None   ,weights=None,
                 ops=None 
                 ):
        """ Constructor for Transformation.
        
        This is designed to have a source and a target.  The dataframes
        will be the same size as the number of elements on the target,
        but will reference elements on the source.

        N is the number of elements in the target mesh.
        M is the number of elements to use in the source mesh.

        Args:
            source (:obj:`~steam.mesh.Mesh()`): Source mesh
            target (:obj:`~steam.mesh.Mesh()`): Target mesh
            elements (:obj:`matrix`): NxM array for source elements to use.
            weights  (:obj:`matrix`): NxM array for corresponding weights.
            ops      (:obj:`list`)  : list of transformation operations performed

        """
        self.source         = source      # This is the source mesh
        self.target         = target      # This is the target mesh
        self.elements       = elements    # Source elements for each Target element
        self.weights        = weights     # weights for each source element.

        self.ops            = ops         # operations (mirror/rotate) to perform on vectors

        self.source_nel     = None
        self.target_nel     = None

        if (source is None or target is None):
            return

        ### source.xyz_el and target.xyz_el need to be defined
        #! Initialize some parameters if I actually passed things in
        addl_error = "\nUnable to get element centroids from {:s} grid!\n"
        try:
            # Number of elements in source
            self.source_nel = self.source.xyz_el.shape[0]   
        except AttributeError as e:
            self.source.get_xyz_el()
            self.source_nel = self.source.xyz_el.shape[0]
            #raise Exception( addl_error.format( 'source' ) )
            #from sys import exc_info
            #raise AttributeError( str(e) + addl_error.format( 'source' ) 
            #                      ).with_traceback( exc_info()[2])

        try:
            # Number of elements in target
            self.target_nel = self.target.xyz_el.shape[0]
        except AttributeError:
            self.target.get_xyz_el()
            self.target_nel = self.target.xyz_el.shape[0]
            #raise Exception( addl_error.format( 'target' ) )
            #from sys import exc_info
            #raise AttributeError( str(e) + addl_error.format( 'target' ) 
            #                      ).with_traceback( exc_info()[2])

        
        #! Make default arrays.  Index is target element number.
        if (self.weights is None and self.elements is None):
            self.weights  = np.zeros([self.target_nel,1])
            self.elements = np.ones( [self.target_nel,1],dtype=np.int_) * -1


    def __str__(self):
        """ Return the print string."""
        string = "Transformation Object:\n"
        return string

    def __repr(self):
        return("STEAM Transformation Object")


    ### ## HDF5 Methods

    def write_hdf5(self,hdf5,root="/",options=None):
        """ Simple routine to write this transformation to disk.

        By default the meshes are written with the transform.
        The 'ioMesh' option can be used to not do this.
        
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write.
                Default is to write to root location "/".
            options (:obj:`dict`, optional): options to vary behavior.

        Current options are:
            ioMesh  = True/False : Should the meshes be written? [T]
        """

        try:
            ioMesh = options['ioMesh']
        except:
            ioMesh = True

        #! I'm just going to pickle or attribute everything
        
        if (ioMesh):
            self.source.write_hdf5(hdf5,root+"/source")
            self.target.write_hdf5(hdf5,root+"/target")

        #!  Make objects that will store attributes
        info = pd.Series(np.zeros(1))
        iroot = root+'/info'
        hdf5.put(iroot,info)
        hdf5.get_storer(iroot).attrs.source_nel      = self.source_nel 
        hdf5.get_storer(iroot).attrs.target_nel      = self.target_nel 

        #!  Pickled things are second since I need to make the groups
        def write_pnode(hdf5,name,obj):
            steam.util.write_hdf5_pnode(
                                   hdf5,
                                   root+"/"+name,
                                   obj,
                                   compress=True
                                   )

        ewo = (self.elements,self.weights,self.ops)
        write_pnode(hdf5,'ewo',ewo)

        return

    def read_hdf5(self,hdf5,root="/",options=None):
        """ Simple routine to read all necessary items from disk.
 
        Source and target meshes are saved by default.  
        'ioMesh' option overwrites the default behavior.
       
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
                Default is to write to root location "/".
            options (:obj:`dict`, optional): options to vary behavior.

        Current options are:
            ioMesh  = True/False : Should the meshes be read? [T]
        """

        try:
            ioMesh = options['ioMesh']
        except:
            ioMesh = True

        if (ioMesh):
            self.target = steam.mesh.Mesh()
            self.target.read_hdf5(hdf5,root+"/target")
            self.source = steam.mesh.Mesh()
            self.source.read_hdf5(hdf5,root+"/source")

        #! I'm just going read the pickles or attributes
        ### Read objects that were stored attributes
        iroot = root+'/info'
        self.source_nel  = hdf5.get_storer(iroot).attrs.source_nel 
        self.target_nel  = hdf5.get_storer(iroot).attrs.target_nel 

        # Read all of the pickles
        def read_pnode(hdf5,name):
            obj = steam.util.read_hdf5_pnode(
                               hdf5,
                               root+"/"+name,
                               compress=True
                               )
            return obj

        (self.elements,self.weights,self.ops) = read_pnode(hdf5,'ewo')

        return


    #! ### Methods
    def set_elements(self,ids,elements):
        """ Set the elements for the transformation.
        
        Size has to conform to the target mesh.
        N is the number of elements in the target mesh.
        M is the number of elements to use in the source mesh.

        Args:
            ids      (:obj:`list`): element indicies (length N)
            elements (:obj:`matrix`): NxM array.
        """

        #! Check dimensions.  Need at least M columns in self.elements
        try:
            new_col = elements.shape[1]
        except IndexError:
            new_col = 1

        old_col = self.elements.shape[1]
        new_col = new_col - old_col
        if (new_col > 0):
            logger.info(" Expanding elements matrix {} columns".format(new_col))
            temp = np.zeros([self.elements.shape[0],new_col],dtype=np.int_)
            self.elements = np.hstack([self.elements,temp])
        if (new_col < 0):
            # Need to update the one I passed in
            new_col = abs(new_col)
            temp = np.zeros([elements.shape[0],new_col])
            elements = np.hstack([elements,temp])

        self.elements[ids] = elements

        #for (i,ind) in enumerate(ids):
        #    if (ind >=self.elements.shape[0]):
        #        raise Exception("Index outside of range in set_elements!")
        #    # This is a sanity check
        #    #if (max(self.elements[ind,:]) > 0):
        #    #    print("Overwriting weights!")

        #    self.elements[ind] = elements[i]

    def set_weights(self,ids,weights):
        """ Set the weights for the transformation.
        
        Size has to conform to the target mesh.
        N is the number of elements in the target mesh.
        M is the number of elements to use in the source mesh.

        Args:
            ids      (:obj:`list`): element indicies (length N)
            weights (:obj:`matrix`): NxM array.
        """

        #! Check dimensions.  Need at least M columns in self.weights
        try:
            new_col = weights.shape[1]
        except IndexError:
            new_col = 1
        old_col = self.weights.shape[1]
        new_col = new_col - old_col
        if (new_col > 0):
            logger.info(" Expanding weights matrix {} columns".format(new_col))
            temp = np.zeros([self.weights.shape[0],new_col])
            self.weights = np.hstack([self.weights,temp])
        if (new_col < 0):
            # Need to update the one I passed in
            new_col = abs(new_col)
            temp = np.zeros([weights.shape[0],new_col])
            weights = np.hstack([weights,temp])
            

        self.weights[ids] = weights

        #for (i,ind) in enumerate(ids):
        #    if (ind >=self.weights.shape[0]):
        #        raise Exception("Index outside of range in set_weights!")
        #    # This is a sanity check
        #    #if (max(self.weights[ind,:]) > 0):
        #    #    print("Overwriting weights!")

        #    self.weights[ind] = weights[i]

    def apply(self,source):
        """ Apply transform to source solution data.

        Use the weights and elements to return a transformed
        solution based off of the source solution passed in.

        Args:
           source (:obj:`~steam.solution.Solution()`): Element-based source solution

        Returns:
           (:obj:`~steam.solution.Solution()`): Element-based target solution
        """


        #! Check to make sure that all elements have been initialized or it will fail.
        if np.min(self.elements) < 0:
            raise Exception("Transform has not been initialized for all elements in target grid!")


        #!  AMS 8/18/17:
        #!    I took this out since this requires that the solution has a mesh.
        #!    that's how these work to check to see if it point or element.  The
        #!    next check is really the one that catches if something is 'wrong'.
        #!if (source.point or source.element is not True):
        #!    raise Exception("Source soln needs to be element based!")

        #! First check to make sure that sizes make sense:
        #! - AMS 9/5/17: I don't think that this is right.  We want to make sure
        #!   that the max value in the elements array is bounded by data shape,
        #!   Not that it has the same number of elemebts as the new soln...
        #if (source.data.shape[0] != self.elements.shape[0]):
        #    print("source.data.shape[0] != self.elements.shape[0]",
        #          source.data.shape[0],
        #          self.elements.shape[0]
        #          )
        #    raise Exception("Source soln does not match source mesh shape!")

        #! Things TODO:
        # print("Add Error check here for NxM size!")
        # print("Add Error check here for vectors/operations!!")

        # OK, so in the case where the indexes are not contiguous on the source data
        # (for instance, if we only include a portion of the geometry).  Then we want
        # to fill the blank indicies with nan in order for numpy's indexing to match.
        # So reindex to the size of the original source that came in.
        # Specific example:
        # - I want to interpolate to a point cloud.  The source data only includes the
        #   elements that are nearest to the point cloud bounds, not all of the source
        #   geometry since most of those points are useless.  The incoming soln includes
        #   a sparce element set.  if I just did a source.data.values, then it would be
        #   an array of len(num_required_elements).  But the .elements list of the
        #   transform assumes the array will be the size of the original grid.

        expanded_df = source.data.reindex(range(self.source_nel))

        source_mat = expanded_df.values


        # PRE 10/17/17 (slow!)
        #for (i,elem) in enumerate(self.elements):
        #   #print(i,elem)
        #    for (j,e) in enumerate(elem):
        #       #print(i,j,e)
        #        target_mat[i]    += source_mat[e] * self.weights[i][j]
        #       #print('s mat ',source_mat[e])
        #       #print('weight',self.weights[i][j])
        #       #print('t mat ',target_mat[i])
        #   #print('t mat ',target_mat[i])

        # POST 10/17/17 (fast!)
        #  - source_mat[elements].shape = [self.target_nel,self.nweights,nq]
        #  - weights.shape = [self.target_nel,self.nweights]
        #  so we want to multiply the center axis and sum the first axis for each of the
        #  third axis
        target_mat = np.einsum('...jk,...j',source_mat[self.elements],self.weights)

        #! Replace the data
        vars   = source.data.columns.tolist()
        new_df = pd.DataFrame(
                              target_mat,
                              columns=vars
                              )

        soln = steam.solution.Solution();
        soln.data = new_df
        soln.store_at = 'ELEMENTS'
        soln.mesh = self.target
        return soln

    @steam.util.timer(logger)
    def inverse_distance_quad(self,k=100,p=2,source_comp=['all'],target_comp=['all']):
        """ Create inverse distance transform between source/target.

        Use an inverse distance method to make a transform.  Uses the
        nearest k points and distance is raised to the p power. This
        routine will ensure that at least one point is from each quadrant,
        so k may be increased to 8 if there is a point in all quadrants.

        It will find the closest in each quadrant - so the maximum number of
        points used is 8 regardless of how high the value of k is that is passed.

        weight = 1/distance^p

        This will populate the transformation matrix in this object.

        Args:
            k (int): Number of nearest points to use when bucketing quadrants.
            p (int): Power to raise distance to.
            source_comp (list): list of components to use in source; Default 'all'
            target_comp (list): list of components to populate in target; Default 'all'

        Return:
            (:obj:`list`): Distance to nearest point for each target point.

        """

        from time import time
        start = time()
        #! Get the elements from the source/target grid
        source_els = self.source.get_comp(source_comp)

        target_els = self.target.get_comp(target_comp)

        try:
            source_index  = np.array(self.source.xyz_el.iloc[source_els].index)
            source_points = self.source.xyz_el.iloc[source_els].values
        except:
            raise Exception("Was not able to get element centroids from Source grid!")

        try:
            target_index  = np.array(self.target.xyz_el.iloc[target_els].index)
            target_points = np.array(self.target.xyz_el.iloc[target_els].values)
        except:
            raise Exception("Was not able to get element centroids from Source grid!")

        #! If there are fewer source points than requested, then lower the k value
        k = min(k,source_points.shape[0])

        logger.debug(" - Startup took {:12f}s".format(time()-start))
        start = time()

        #! Build the KDTree
        kd_tree = cKDTree(source_points)

        logger.debug(" - KD took {:12f}s".format(time()-start))
        start = time()

        #! Get the nearest k points to each target point
        ele_near = kd_tree.query(target_points,k=k)
        # ele_near[0] is distance matches
        # ele_near[0][0] are the distances to the K closests elements to element 0
        # ele_near[1][0] are the K closests elements to element 0
        (dist,index) = ele_near

        if k == 1:
            #! Convert ele_near indexing into global indexing:
            elements = source_index[np.array(index)]
    
            #! Need to make all of these 1D arrays 2D
            weights  = np.array([[1] for i in elements])
            elements = np.array([[i] for i in elements])

            #! Transform object needs elements and weights
            self.set_weights( target_index,weights)
            self.set_elements(target_index,elements)

            return dist

        #! Return the distance to the nearest point for each target point
        dist_return = dist[:,0]

        #! This is the distance from the target_point to the neighbor
        #vectors = source_points[index][:] - target_points[:,np.newaxis,:]
        vectors = source_points[index,:] - target_points[:,np.newaxis,:]
        # Shape (num_tp,num_k,3)

        logger.debug(" - Lookup took {:12f}s".format(time()-start))
        start = time()

        #! Now get subsets for each quadrant
        # - odd trick to get my < and > make sense:
        #print (i,i%2 == 0,i%4<2,i<4)
        # 0 True  True  True
        # 1 False True  True
        # 2 True  False True
        # 3 False False True
        # 4 True  True  False
        # 5 False True  False
        # 6 True  False False
        # 7 False False False
        quad_dist  = []
        quad_index = []
        for i in range(8):
            if (i%2 == 0):
                b1 = vectors[:,:,0] >= 0
            else:
                b1 = vectors[:,:,0] <  0

            if (i%4 <  2):
                b2 = vectors[:,:,1] >= 0
            else:
                b2 = vectors[:,:,1] <  0

            if (i   <  4):
                b3 = vectors[:,:,2] >= 0
            else:
                b3 = vectors[:,:,2] <  0

            q1     = np.array([b1,b2,b3])
            qindex = np.select([q1.all(axis=0)],[index],default=-1)
            qdist  = np.select([q1.all(axis=0)],[dist ],default= 0)
            qfirst = np.argmax(qindex>=0,axis=1)
            qi     = qindex[np.arange(qindex.shape[0]),qfirst]
            qd     = qdist[ np.arange(qindex.shape[0]),qfirst]
            quad_dist.append(qd)
            quad_index.append(qi)

        new_dist = np.einsum('ij->ji',np.array(quad_dist           ))
        new_index= np.einsum('ij->ji',np.array(quad_index,dtype=int))

        new_dist[ new_index < 0] = np.inf
        new_index[new_index < 0] = 0.0

        logger.debug(" - Quads took {:12f}s".format(time()-start))
        start = time()

        dist  = new_dist
        index = new_index

        #! I want to divide by the array; if the array has a zero value, then
        #! set the weight equal to 1.0
        with np.errstate(divide='ignore'):
            weights = 1.0 / dist ** p

        #! Where do we have the exactly correct point?
        #! At this point, dist has 'inf' values where nothing was found in that quad
        #! If it has any zero values, then we want to inf everything that is non-zero
        #! in that row.  Then the weight should be set to 1.0 for the remaining point.
        exact_array = np.where(dist == 0)
            # Zero the rows that contain an exact match
        weights[exact_array[0],:] = 0.0
            # Set the match locations back to 1.0
        weights[exact_array] = 1.0

        #weights[np.invert(np.isfinite(weights))] = 1.0  # THIS IS WRONG

        #! Sum the weights for each element
        row_sums = np.sum(weights,axis=1)

        #! Broadcast a divide by the sum:
        weights = np.einsum('ij,i->ij',weights,1/row_sums)

        #! Convert ele_near indexing into global indexing:
        elements = source_index[index]

        #! Transform object needs elements and weights
        self.set_weights( target_index,weights)
        self.set_elements(target_index,elements)

        return dist_return
        

    @steam.util.timer(logger)
    def inverse_distance(self,k=5,p=2,source_comp=['all'],target_comp=['all']):
        """ Create inverse distance transform between source/target.

        Use an inverse distance method to make a transform.  Uses the
        nearest k points and distance is raised to the p power.

        weight = 1/distance^p

        This will populate the transformation matrix in this object.

        Args:
            k (int): Number of nearest points to use.
            p (int): Power to raise distance to.
            source_comp (list): list of components to use in source; Default 'all'
            target_comp (list): list of components to populate in target; Default 'all'

        Return:
            (:obj:`list`): Distance to nearest point for each target point.

        """

        #! Get the elements from the source/target grid
        source_els = self.source.get_comp(source_comp)

        target_els = self.target.get_comp(target_comp)

        try:
            source_index  = np.array(self.source.xyz_el.iloc[source_els].index)
            source_points = self.source.xyz_el.iloc[source_els].values
        except:
            raise Exception("Was not able to get element centroids from Source grid!")

        try:
            target_index  = np.array(self.target.xyz_el.iloc[target_els].index)
            target_points = self.target.xyz_el.iloc[target_els].values
        except:
            raise Exception("Was not able to get element centroids from Source grid!")

        #! If there are fewer source points than requested, then lower the k value
        k = min(k,source_points.shape[0])

        #! Build the KDTree
        kd_tree = cKDTree(source_points)

        #! Get the nearest k points to each target point
        ele_near = kd_tree.query(target_points,k=k)
        # ele_near[0] is distance matches
        # ele_near[0][0] are the distances to the K closests elements to element 0
        # ele_near[1][0] are the K closests elements to element 0
        (dist,index) = ele_near

        if k == 1:
            #! Convert ele_near indexing into global indexing:
            elements = source_index[np.array(index)]
    
            #! Need to make all of these 1D arrays 2D
            weights  = np.array([[1] for i in elements])
            elements = np.array([[i] for i in elements])

            #! Transform object needs elements and weights
            self.set_weights( target_index,weights)
            self.set_elements(target_index,elements)

            return dist

        #! I want to divide by the array; if the array has a zero value, then
        #! set the weight equal to 1.0
        with np.errstate(divide='ignore'):
            weights = 1.0 / dist ** p

        #! If it is not finite, set it to 1.0
        #! Where do we have the exactly correct point
        exact = np.any(np.invert(np.isfinite(weights)),1)
        weights[exact,:] = 0.0
        weights[exact,0] = 1.0
        #weights[np.invert(np.isfinite(weights))] = 1.0  # THIS IS WRONG

        #! Sum the weights for each element
        row_sums = np.sum(weights,axis=1)
        #! Broadcast a divide by the sum:
        weights = np.einsum('ij,i->ij',weights,1/row_sums)

        #! Convert ele_near indexing into global indexing:
        elements = source_index[np.array(index)]

        #! Transform object needs elements and weights
        self.set_weights( target_index,weights)
        self.set_elements(target_index,elements)

        #! Return the distance to the nearest point for each target point
        return dist[:,0]
        

class IntSpace():
    """Class definition for STEAM Interpolation Space

    An interpolation space is an object used for flightspace
    interpolation.  It is designed to take an entire or subset of
    a database or np.ndarray.  Data sets can have multiple interpolation 
    spaces, but the reverse is not true.
    """

    def __init__(self,indep=None,name="IntSpace",
                 subvars=None, subvals=None, subtypes=None,
                 subspace=False):
        """Constructor for Interpolation Space.

        This just defines a name and list of independant variables as well as parameters
        for any subspaces.

        There are currently two types of interpolation space:
            * `tri`: uses scipy and qhull to create a simplex-filled Delaunay space
            * `1d` : uses scipy to create a 1d interpolation space


        Args:
            indep       : Independant parameters
            name        : Name of the space
            subvars     : variables to cast subspaces into
            subvals     : indicies in subvars to use for subspaces
            subtypes    : type of interpolation to use for subspaces [Not Implimented]
            subspace    : Is this a subspace?  Internal class use, only.

        """

        self.out_index  = None          # Index conversion to output (in case I was sparce in data set)

        self.name      = name           # Name of the IntSpace

        self.scale     = dict()         # Scale factor for normalization
        self.opts      = dict()         # Mathematical operations to perform on indep params
        self.opt_names = dict()         # String representation of the options

        if indep is None:
            indep = []
        self.indep     = indep          # independent parameters

        ### https://github.com/scipy/scipy/blob/7252c1f3bc0643daa81ff27b2b615181906cfc99/scipy/spatial/qhull.pyx#L1693
        self.space       = None           # Delaunay Triangulation

        self.type      = '?'

        ### Subspaces are Interpolation Spaces that this space interpolates between.  Originally
        ### implimented for ND+M interpolation.
        self.subspaces = None
        self.subvars   = None
        self.subvals   = None
        self.subopts   = None

        self.is_subspace = subspace

        if subvars is None:
            subvars = []

        # If subvars were passed in, then I have subspaces
        if len(subvars) > 0:
            self.type      = '1d'
            self.subspaces = dict()
            self.subvars   = list(subvars)
            self.subopts           = dict()
            self.subopts['vals']   = dict()
            self.subopts['type']   = dict()
            # Ensure that all of the subvars are in the idep list
            for var in subvars:
                if var not in indep:
                    raise KeyError("Subvariable '{}' not an independant parameter.".format(var))
                # Default subvals and subtypes
                self.subopts['vals'][var] = None
                self.subopts['type'][var] = self.type

            # If subvals is specified for any of the keys, overwrite the default
            if subvals is not None:
                for (key,val) in subvals.items():
                    if val is not None:
                        self.subopts['vals'][key] = np.array(val)
            # If subtype is specified, overwrite the default
            if subtypes is not None:
                for (key,val) in subtypes.items():
                    self.subopts['type'][key] = val



    def __str__(self):
        """Return the print string of the IntSpace object."""

        string  = "\nInterpolation Space Object: {}\n".format(self.name)

        string += " Independent Variables  : "+" ".join(self.indep)+"\n"

        string += " Type  : {}\n".format(self.type)
        
        if (self.type == 'tri'):
            if (self.space  is not None):
                string += " Triangulation : Yes"+"\n"
            else:
                string += " Triangulation : None"+"\n"
        
        if (self.type == '1d'):
            if (self.space  is not None):
                string += " 1D Space : Yes"+"\n"
            else:
                string += " 1D Space : None"+"\n"

        string += " Variable operations:\n"
        for (i,key) in enumerate(self.indep):
            if key in self.opt_names:
                val = self.opt_names[key]
                string+="   {}: {}\n".format(key,val)

        string += " Scaling factors:\n"
        for (i,key) in enumerate(self.indep):
            string +="   {}: {:5g}\n".format(key,self.scale[key])

        
        string += " Subspaces :"
        if self.subspaces  is None:
            string +="None\n"
        else:
            string += " ".join([ "{:5g}".format(v) for v in sorted(self.subspaces.keys())])+"\n"
            string+="   Subvals :"+" ".join([ "{:5g}".format(v) for v in sorted(
                                                    self.subvals / self.scale[self.subvars[0]])]
                                                                                )+"\n"

            string+="   Subvars    : {}\n".format(str(self.subvars))
            string+="   Non-subvars: {}\n".format(str(list(
                                            set(self.indep).difference(self.subvars)
                                                  )))
            string+="   Summary: \n"
            string+=self.subspace_string(indent=4)


        return string
 
    def subspace_string(self,indent=0):
        """ Utility to roll up and format subspace vals for all of my subspaces.

        If I don't have subspaces, then I return the dimensionality of my space.

        This is recursive.

        """
        string = ""
    
        space = "".join([" " for i in range(indent)])

        if self.subspaces is None:
            return string

        subvar   = self.subvars[0]
        subscale = self.scale[subvar]
        ## If I am the last subvar, then print a summary
        if len(self.subvars) == 1:
                vals = " ".join(["{:.5g}".format(v / subscale) for v in self.subvals])
                string += space+"{} : {}\n".format(subvar,vals)

        ## Otherwise, print the hierarchy
        else:
            for val in sorted(self.subvals):
                string += space+"{} : {:.5g}\n".format(subvar,val / subscale)
                string += self.subspaces[val].subspace_string(indent+1)

        return string



    def __repr__(self):
        return "STEAM IntSpace Object"

    #! This is just something to give the lambda scope.
    #! Otherwise, variables can clobber each other
    def create_lambda_opt(self,opt,key,verbose=True):
        """ Creates a lambda function operation.

        This keys off of text operations because we need to be able to
        save these with pickle.  Default pickle cannot save lambda functions
        so we save the strings.  If we want to save the lambda functions, then
        we need to 'import dill as pickle', but that would require a new
        dependancy (dill).

        Args:
            opt (:obj:`str`): Option to perform
            key (:obj:`str`): Variable to operate on.
            verbose (:obj:`bool`): Output to screen? [T]

        Returns:
            :obj:`lambda`: Lambda function
        """

        import math

        if opt == "log10":
            f  = lambda x: math.log10(x)
            if verbose:
                print("    - {} -> {}".format(key,'log_10({})'.format(key)))
        elif opt == "loge":
            f  = lambda x: math.log(x)
            if verbose:
                print("    - {} -> {}".format(key,'log_e({})'.format(key)))
        elif '**' in opt:
            pow = float(opt.replace('**',''))
            f  = lambda x: x**pow
            if verbose:
                print("    - {} -> {}".format(key,'{}**{}'.format(key,pow)))
        else:
            raise TypeError("Illegal option for opt '{}'!".format(opt))

        return f


    ### ## HDF5 Methods

    def write_hdf5(self,hdf5,root="/"):
        """ Simple routine to write all of this interpolation space to disk.
        
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write.
                Default is to write to root location "/".
        """

        #! I'm just going to pickle or attribute everything

        #!  Make objects that will store attributes
        info = pd.Series(np.zeros(1))
        iroot = root+'/info'
        hdf5.put(iroot,info)
        hdf5.get_storer(iroot).attrs.name        = self.name
        hdf5.get_storer(iroot).attrs.indep       = self.indep
        hdf5.get_storer(iroot).attrs.scale       = self.scale
        hdf5.get_storer(iroot).attrs.opt_names   = self.opt_names
        hdf5.get_storer(iroot).attrs.type        = self.type
        hdf5.get_storer(iroot).attrs.is_subspace =  self.is_subspace

        #!  Pickled things are second since I need to make the groups
        steam.util.write_hdf5_pnode(
                               hdf5,
                               root+"/space",
                               self.space,
                               compress=True
                               )

        steam.util.write_hdf5_pnode(
                               hdf5,
                               root+"/out_index",
                               self.out_index
                               )
        #! Subspace info
        if self.subspaces is None:
            subspaces = []
        else:
            subspaces = [*self.subspaces.keys()]  # List of keys
        steam.util.write_hdf5_pnode(
                               hdf5,
                               root+"/subspace_info",
                                (subspaces,
                                 self.subvars,
                                 self.subvals,
                                 self.subopts
                                )
                               )
        for (s,sub) in enumerate(subspaces):
            self.subspaces[sub].write_hdf5(hdf5,
                                           root+"/subspace_{:d}".format(s)
                                           )

        return

    def read_hdf5(self,hdf5,root="/"):
        """ Simple routine to read all necessary items from disk.
        
        Args:
            hdf5 (`~pandas.io.pytables.HDFStore`): HDFStore object.
            root (:obj:`str`,optional): Location in HDF5 to write database.
                Default is to write to root location "/".
        """

        #! I'm just going to pickle or attribute everything
        
        self.space    = steam.util.read_hdf5_pnode(
                               hdf5,
                               root+"/space",
                               compress=True
                               )

        self.out_index = steam.util.read_hdf5_pnode(
                               hdf5,
                               root+"/out_index",
                               )
        #! Subspace info
        (subspaces,
         self.subvars,
         self.subvals,
         self.subopts) = steam.util.read_hdf5_pnode(
                               hdf5,
                               root+"/subspace_info",
                               )
        if len(subspaces) > 0:
            self.subspaces = dict()
            for (s,sub) in enumerate(subspaces):
                self.subspaces[sub] = steam.interpolate.IntSpace()
                self.subspaces[sub].read_hdf5(hdf5,
                                               root+"/subspace_{:d}".format(s)
                                               )

         ### Make objects that will store attributes
        iroot = root+'/info'
        self.name  = hdf5.get_storer(iroot).attrs.name
        self.indep = hdf5.get_storer(iroot).attrs.indep

        try: 
            self.is_subspace = hdf5.get_storer(iroot).attrs.is_subspace
        except:
            self.is_subspace = False

        try:
            self.scale = hdf5.get_storer(iroot).attrs.scale
        except:
            self.scale = dict()

        try:
            self.type  = hdf5.get_storer(iroot).attrs.type
        except:
            self.type  = 'tri'

        #! Load/create the lambda functions
        try:
            self.opt_names = hdf5.get_storer(iroot).attrs.opt_names
        except:
            self.opt_names = dict()

        for (key,opt) in self.opt_names.items():
            f = self.create_lambda_opt(opt,key,verbose=False)
            self.opts[key]      = f



        return

    ### ## Intepolation Methods


    def build_space(self,points,out_index=None,scale=None,opts=None):
        """Subroutine to build the interpolation space and any relevant subspaces.

        Points is an ndarray of points.  The columns need to line-up with the
        independant parameters that were passed in at creation.

        If the 'points' array is a subset of your larger array, you can pass the
        indicies (in the larger array) for the subset as 'out_index'.  When called
        to return interpolation weights and stencil, the stencil will be upated to
        conform to the values in out_index.  If not provided, then the row index
        in 'points' will be returned.

        scale is a dict where the keys desribe what to multiply axes by for scaling.
        If scale is the string 'auto' then do auto-scaling to all parameters.  Scaling
        is applied after any manipulations specified in 'opts'.
        
        Examples:
        * scale='auto' - Auto-scale each paramters to be unit range
        * scale={'p1':'auto','p2':1.0} - Auto-scale parameter 'p1', don't scale 'p2'

        opts is a dict where manipulations of the independant paramters can be performed
        to change the specifics for the linear interpolation. Possible opts are the following:
            * "loge"   - Scale by log_e
            * "log10"  - Scale by log_10
            * "\*\*X"    - Raise to the X power
            * lambda function - provide your own operation [TODO!]

        Examples:
          * opts={'p1':"\*\*3",'p3':"log10"} - Perform interpolation in 3-space: (p1^3, p2, log10(p3))

        Args:
            points (:obj:`np.ndarray`): Point locations in space.
            out_index (:obj:`np.ndarray`): index mapping for points data.
            scale (:obj:`dict`): Optional scaling parameters to perform on the indep params.
            opts (:obj:`dict`): Optional operations to perform on the indep params.

        """

        ### Check the length of the independent/dependent variables
        if self.indep is None or len(self.indep) == 0:
            raise Exception("ERROR: Indepenent variables are not defined!")

        ### Specify index mapping
        if out_index is not None:
            self.out_index = np.array(out_index,dtype="int")
        else:
            self.out_index = np.arange(points.shape[0])
        
        #! Promote 1D arrays to 2D columns
        try:
            ndim = points.shape[1]
        except IndexError:
            points = np.array([points]).T
            ndim = 1

        #! Copy the points to a new array since I'm going to manipulate it
        spoints = np.array(points,dtype='float')
        
        #! Subspaces don't have to normalize anything because it is passed in normalized
        if self.is_subspace:
            self.scale = scale
            self.opts  = opts

        else:
            #! Do we need to operate on the point data?
            if opts is not None:
                print(" - Variable Operations:")
                for key in self.indep:
                    if key in opts:
                        opt = opts[key]
                        f = self.create_lambda_opt(opt,key)
                        self.opts[key]      = f
                        self.opt_names[key] = opt
                        
            #! Evaluate the math operations
            for (i,key) in enumerate(self.indep):
                if key in self.opts:
                    #print(key,2,self.opts[key](2))
                    spoints[:,i] = np.array([self.opts[key](x) for x in spoints[:,i]])


            #! Initialize the scaling parameters
            auto_scale = dict()
            no_scale   = dict()
            for (i,key) in enumerate(self.indep):
                var_range       = max(spoints[:,i]) - min(spoints[:,i])
                auto_scale[key] = 1.0 / var_range
                no_scale[key]   = 1.0

            if scale == 'auto' :
                self.scale = auto_scale
            elif type(scale) is type(dict()):
                for (i,key) in enumerate(self.indep):
                    if key in scale:
                        if (scale[key] == "auto"):
                            self.scale[key] = auto_scale[key]
                        else:
                            self.scale[key] = scale[key]
                    else:
                        self.scale[key] = no_scale[key]
            else:
                self.scale = no_scale
                

            #! Evaluate the scaling parameters
            print(" - Scaling factors:")
            for (i,key) in enumerate(self.indep):
                if key in self.scale:
                    spoints[:,i] *= self.scale[key]
                    print("    - {}: {}".format(key,self.scale[key]))



        # If I don't have any subspaces, then make the Delaunay triangulation
        if self.subspaces is None:

            if ndim == 1:
                self.type = "1d"
                # Pass as transpose to make it 1D
                self.build_linear(spoints.T[0])
            else:
                self.type = "tri"
                self.build_delaunay(spoints)

        # If I have subspaces, then I'm going to make each of them now
        else:

            subvar  = self.subvars[0]
            subvari = self.indep.index(subvar)

            # Determine the number of values and therefore subspaces
            self.subvals = np.unique(spoints[:,subvari])
 
            # If specific subspace values were passed in, then use them.  However,
            # only include the subspace values if they exists in this space.
            if self.subopts['vals'][subvar] is not None:
                self.subvals = np.intersect1d(self.subvals,
                # If I passed in explicit values in __init__, I need to scale them to be consistant
                                        self.subopts['vals'][subvar] * self.scale[subvar]
                                        )

            # Just in case they weren't sorted
            self.subvals = np.sort(self.subvals)

            self.subspaces = dict()
            for (i,subval) in enumerate(self.subvals):
                dsubval = subval / self.scale[subvar]
                logger.debug("making subspace for {} == {}".format(subvar,dsubval))
                
                subpoints = np.delete(spoints[spoints[:,subvari] == subval],subvari,1)
                subindex  = self.out_index[spoints[:,subvari] == subval]


                # Subspaces are not going to have the subvar as a variableyy
                subindep = list(self.indep)
                subindep.remove(subvar)

                # Subspaces any remaining subvars
                subsubs  = self.subvars[1:]

                name = "{}.{}_{}".format(self.name,subvar,dsubval)
                self.subspaces[subval] = IntSpace(
                                        indep    = subindep,
                                        name     = name,
                                        subvars  = subsubs,
                                        subvals  = self.subopts['vals'],
                                        subtypes = self.subopts['type'],
                                        # I have to declare it a subspace to turn
                                        # scaling off and other preprocessing steps
                                        subspace = True,
                                        )

                #print(subvar,i,subval)
                # Scaling and opts are already done, so just pass the points
                self.subspaces[subval].build_space(
                                        subpoints,
                                        out_index = subindex,
                                        # Pass the scale so that I can translate between
                                        # dimensional subvals and the non-d points
                                        scale     = self.scale,
                                        # Pass in the opts just to have them
                                        opts      = self.opts,
                                        )

        return

    def build_linear(self,points=None):
        """Subroutine to build the :obj:`~scipy.interp1d` object.

        Points is an ndarray of points.   It needs to be one-dimensional.

        Args:
            points (:obj:`np.ndarray`): Point locations in space.

        """

        ### Create the 1D space triangulation
        import scipy.interpolate as interpolate

        ### I need to sort it and the verticies
        order = np.argsort(points)
        self.out_index = self.out_index[order]
        points = points[order]

        self.space = interpolate.interp1d(points,np.arange(len(points)))

        print(" - Built {}-D interpolation space".format(len(self.indep)))
            
        return

    def build_delaunay(self,points=None):
        """Subroutine to build the :obj:`~scipy.Delaunay` triangulation.

        Points is an ndarray of points.  The columns need to line-up with the
        independant parameters that were passed in at creation.

        Args:
            points (:obj:`np.ndarray`): Point locations in space.

        """

        ### Create the Delaunay triangulation
        self.space = Delaunay(points)

        logger.info(" - Built {}-D Delaunay interpolation space".format(len(self.indep)))

        if (self.space.coplanar.size > 0):
            print(" --> WARNING: Not all verticies were used!")
            print(" --> {} unused.  Try running Qhull with QJ.".format(self.space.coplanar.shape[0]))
            
        return



    def get_weights(self,  cases):
        """Method that returns the indicies and weights to use for interpolation.

        Linear interpolation weights are returned for each point provided.
        
        Args: 
            cases (:obj:`DataFrame`, :obj:`list`, or :obj:`~numpy.ndarray`): Indep parameters to return weights for.

        Returns:
            :obj:`list`: List of vertex indicies for each interpolation point
            :obj:`list`: List of vertex weights for each interpolation point
        """

        ## We can pass in a dataframe or a numpy array
        if type(cases) == pd.DataFrame:
            #### Check to make sure that the caselist has the independent 
            #### variables that this space depends on.
            #if not steam.util.vars_in_list(
            #        self.indep,
            #        cases.keys()
            #        ):
            #    print("Case: ",cases.keys())
            #    print("IntSpace: ",self.indep)
            #    raise Exception("Not all database vars included in cases!")

            ### Make a point that has the variables lined up correctly
            points = cases[self.indep].values

        else:
            # Assume it is a numpy array with the correct variables and order
            points = np.array(cases,dtype='float')
        
        #! Promote 1D arrays to 2D columns
        try:
            ndim = points.shape[1]
        except IndexError:
            points = np.array([points],dtype='float').T

        try:
            nrow = points.shape[0]
            ncol = points.shape[1]
        except:
            logging.error("Input does not look like a DataFrame or numpy ndarray")
            raise
        if ncol != len(self.indep):
            raise IndexError("Number of columns is not the same as the number of independant params")
        if nrow <= 0:
            raise IndexError("I require at least one row/case to get weights for")

        #! Done input checking

        if not self.is_subspace:
            #! Evaluate the math operations
            for (i,key) in enumerate(self.indep):
                if key in self.opts:
                    points[:,i] = np.array([self.opts[key](x) for x in points[:,i]])

            #! Evaluate the scaling parameters
            scale_arr = np.array([self.scale[var] for var in self.indep])
            points *= scale_arr
            #if self.scale is not None:
            #    for (i,key) in enumerate(self.indep):
            #        if key in self.scale:
            #            points[:,i] *= self.scale[key]



        # If I don't have any subspaces, then return the Delaunay weights
        if self.subspaces is None:

            if self.type == 'tri':
                (verts,weights) = self.get_weights_delaunay(points)
            elif self.type == '1d':
                (verts,weights) = self.get_weights_linear(points.T[0])
            else:
                raise ValueError("Unknown IntSpace type: {}".format(self.type))

        # If I have subspaces, then I'm going to make each of them now
        else:

            (verts,weights) = self.get_weights_subspace(points)
        
        return (verts,weights)



    def get_weights_subspace(self,  points):
        """Method that returns the weights using results from the subspaces.

        Linear interpolation weights are returned for each point provided.  The
        point indicies are in the out_index convention.

        Args: 
            points (:obj:`np.array`): Array of points in space to interpolate to.

        Returns:
            :obj:`list`: List of vertex indicies for each interpolation point
            :obj:`list`: List of vertex weights for each interpolation point
        """

        subvar  = self.subvars[0]
        subvari = self.indep.index(subvar)

        weights = []
        verts   = []
        for (i,point) in enumerate(points):

            ### Find the location in the subspace variable and put it
            ### back into normal space
            subval = point[subvari]
            ### Make the point a 2D array (one row) without my subvar
            subpoint  = np.array([np.delete(point,subvari)])
        
            ### Am I inside of the subspace range?
            if subval < self.subvals[0] or subval > self.subvals[-1]:
                vert   = self.out_index[:1]
                weight = np.array([-1.0])
                #raise Exception("Out of bounds ({}): {:2g} [{:.2g} {:.2g}]".format(
                #                                            self.name,
                #                                            subval,
                #                                            self.subvals[0],
                #                                            self.subvals[-1]
                #                                            ))
            else:
                # Find the bracket of cases
                val_loc = self.subvals.searchsorted(subval)
                #print(self.subvals)
                #print(self.subvals <= subval)

                sval_min = self.subvals[max(val_loc - 1,0)]
                sval_max = self.subvals[val_loc]

                #print("{:10s} values: {:2g} [{:.2g} {:.2g}]".format(
                #                    subvar,
                #                    subval,
                #                    sval_min,
                #                    sval_max))

                if subval == sval_max:
                    # Just do the interpolation and return
                    subspace = self.subspaces[sval_max]
                    (vert, weight) = subspace.get_weights(subpoint)
                    vert   = vert[0]
                    weight = weight[0]
                else:
                    # I need to bracket the interpolation
                    subspace = self.subspaces[sval_min]
                    (vert_min, weight_min) = subspace.get_weights(subpoint)

                    subspace = self.subspaces[sval_max]
                    (vert_max, weight_max) = subspace.get_weights(subpoint)

                    max_rat = (subval - sval_min) / (sval_max - sval_min)
                    min_rat = 1.0 - max_rat
                    #print("Subvar weights:     [{:.2g} {:.2g}]".format(
                    #                    min_rat,
                    #                    max_rat))

                    vert   = np.concatenate((vert_min[0],vert_max[0]))
                    weight = np.concatenate((weight_min[0] * min_rat,
                                             weight_max[0] * max_rat))
                    #print(vert)
                    #print(weight)

            weights.append(weight)
            verts.append(vert)

        return (verts,weights)

    def get_weights_linear(self,  points):
        """Method that returns the weights using the 1D interpolation space.

        Linear interpolation weights are returned for each point provided.  The
        point indicies are in the out_index convention.

        Args: 
            points (:obj:`np.array`): Array of points in space to interpolate to.

        Returns:
            :obj:`list`: List of vertex indicies for each interpolation point
            :obj:`list`: List of vertex weights for each interpolation point
        """

        
        weights = []
        verts   = []
        for p in points:

            ### Find the indexes from the interpolation space
            try:
                idx  = self.space(p)
            except ValueError:
                # Out of bounds!
                vert   = [   0]
                weight = [-1.0]

            else:
                ### The index needs to be turned into weights
                ### For example, 1.5 -> 0.5 weight for index 1 and 0.5 weight for index 2
                ###              1.3 -> 0.7 weight for index 1 and 0.3 weight for index 2
                ###              0.3 -> 0.7 weight for index 0 and 0.3 weight for index 1
                np1 = int(np.floor(idx))
                np2 = int(np.ceil(idx))

                if np1 == np2:
                    vert   = [np1]
                    weight = np.array([1.0])
                else:
                    w2 = (idx - np1) / (np2 - np1)
                    w1 = 1.0 - w2
                    vert   = [np1,np2]
                    weight = np.array([ w1, w2])

            weights.append(weight)

            #! Verts need to be put back into the index space of the parent database
            vert = self.out_index[vert]

            verts.append(vert)

        return (verts,weights)

    def get_weights_delaunay(self,  points):
        """Method that returns the weights using Delaunay interpolation space.

        Linear interpolation weights are returned for each point provided.  The
        point indicies are in the out_index convention.


        Args: 
            points (:obj:`np.array`): Array of points in space to interpolate to.

        Returns:
            :obj:`list`: List of vertex indicies for each interpolation point
            :obj:`list`: List of vertex weights for each interpolation point
        """

        ### Find the simplex that this is in
        simps  = self.space.find_simplex(points)

        weights = []
        verts   = []
        for (i,simp) in enumerate(simps):

            point = points[i,:]
            ### Get the verts of the simplex (indicies in self.data)
            vert    = self.space.simplices[simp]

            # T c = x - r
            # where c are the barycentric coordinates
            #       x is the test point
            #       r is the vector distance to the ndim+1 vertex 
            #            (self.space.transform[simplex,ndim])
            #       T^(-1) = self.space.transform[simplex,:ndim] - inverse of T
            # So to get the weights for 1-ndim, use c = T(x-r).  
            # The weight for ndim+1 is from knowing that all of the 
            #   barycentric coordinates must equal one.

            ndim = len(self.indep)
            b    = self.space.transform[simp,:ndim].dot(point - self.space.transform[simp,ndim])
            w    = np.append(b,1.0-sum(b))  # Force the sum to be 1.0

            # Is this the best place for this?
            import math
            if (math.isnan(b[0])):
                w[:] = -1.0

            weights.append(w)

            #! Verts need to be put back into the index space of the parent database
            vert = self.out_index[vert]

            verts.append(vert)

        return (verts,weights)


###############################################################################
#!  ## Functions

def write_tec_tets(intspace, fname='delaunay.dat'):
    """For databases of three dimensions or less, write a *.dat file that 
    contains connectivity and can be plotted in Tecplot.
    This can be useful to investigate how the Delaunay triangulation
    performed.
    
    Args:
        intspace (:obj:`~steam.interpolate.IntSpace`): Interpolation space to be plotted.
        fname (:obj:`str`): Output file name. 
            Default is 'delaunay.dat'.
    """

    assert len(intspace.indep) <= 3, (
        'Interpolation space has > 3 independent variables.  Cannot plot.' )

    if steam.has_dataio:
        import pydata
    
        tec = pydata.data.PyData ()
    
        tec.xyz_vars  = intspace.indep
    
        tec.nodes     = intspace.space.points.transpose ()
        for i, key in enumerate (intspace.indep):
            tec.nodes [i] /= intspace.scale [key]
    
        tec.connectivity = (intspace.space.simplices).tolist ()
    
        tec.components   = [1 for i in range (intspace.space.simplices.shape[0])]
    
        tec.write_file (fname, 'tec')

    else:

        ### Assemble the file header
        header_1 = 'TITLE="None"\nVARIABLES ='
        for var in intspace.indep:
            header_1 += ' "{:s}"'.format(var)
        header_2 = ( 'ZONE T="IntSpace"\nZONETYPE="FETETRAHEDRON"\n' + 
                     'DATAPACKING="BLOCK"\nN={:d} E={:d}'.format(  
                     intspace.tri.points.shape[0], intspace.tri.nsimplex ) )
        head = "{:s}\n{:s}".format( header_1, header_2 )
    
        ### Write to the output file
        #   Write point coordinates
        with open( fname, 'wb' ) as f:
            np.savetxt( f, intspace.tri.points[:,0], header=head, comments='' )
        with open( fname, 'ab' ) as f:
            np.savetxt( f, intspace.tri.points[:,1] )
            np.savetxt( f, intspace.tri.points[:,2] )
    
            ### Save connectivity -- add 1 to every node number
            np.savetxt( f, intspace.tri.simplices + 1 )

def write_vtu_tets(intspace, fname='delaunay.vtu'):
    """For databases of three dimensions or less, write a vtu file that 
    contains connectivity and can be plotted via external tools.
    This can be useful to investigate how the Delaunay triangulation
    performed.
    
    Args:
        intspace (:obj:`~steam.interpolate.IntSpace`): Interpolation space to be plotted.
        fname (:obj:`str`): Output file name. 
            Default is 'delaunay.vtu'.
    """
    ### Modified from :
    ###  https://github.com/cfinch/Shocksolution_Examples/blob/master/Visualization/vtktools.py

    import xml.dom.minidom

### <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
### <UnstructuredGrid>
### <Piece NumberOfPoints="76" NumberOfCells="14">
### <CellData>
### <DataArray type="Float64" Name="r" Format="ascii">
### 0.100000E+001
### ...
### </DataArray>
### <Points>
### <DataArray NumberOfComponents="3" type="Float64" Format="ascii">
###   4.6666666666666661       -1.0000000000000000        1.0000000000000000
### ...
### </DataArray>
### </Points>
### <Cells>
### <DataArray type="Int32" Name="connectivity" Format="ascii">
###    0     1     2     3     4     5     6      8     9    10    11    12    13    14
### ...
### </DataArray>         
### <DataArray type="Int32" Name="offsets" Format="ascii">
###          8          
###           16
### </DataArray>
### <DataArray type="Int32" Name="types" Format="ascii">
###          12
### </DataArray>
### </Cells>
### </Piece>
### </UnstructuredGrid>
### </VTKFile>

    doc = xml.dom.minidom.Document()
    root_element = doc.createElementNS("VTK", "VTKFile")
    root_element.setAttribute("type", "UnstructuredGrid")
    root_element.setAttribute("version", "0.1")
    root_element.setAttribute("byte_order", "LittleEndian")
    doc.appendChild(root_element)

    # Unstructured mesh element
    unstructuredGrid = doc.createElementNS("VTK", "UnstructuredGrid")
    root_element.appendChild(unstructuredGrid)

    # Piece 0 (only one)
    piece = doc.createElementNS("VTK", "Piece")
    piece.setAttribute("NumberOfPoints", str(intspace.space.points.shape[0]))
    piece.setAttribute("NumberOfCells", str(intspace.space.nsimplex))
    unstructuredGrid.appendChild(piece)

    ### Points ####
    points = doc.createElementNS("VTK", "Points")
    piece.appendChild(points)

    # Point location data
    point_coords = doc.createElementNS("VTK", "DataArray")
    point_coords.setAttribute("type", "Float64")
    point_coords.setAttribute("format", "ascii")
    point_coords.setAttribute("NumberOfComponents", "3")
    points.appendChild(point_coords)

    string = "\n"
    for i in range(intspace.space.points.shape[0]):
        point = intspace.space.points[i,:]

        #! Undo the scaling parameters
        if type(intspace.scale) is type(dict()):
            for (j,key) in enumerate(intspace.indep):
                if key in intspace.scale:
                    point[j] /= intspace.scale[key]

        string += " {:.12f} {:.12f} {:.12f}\n".format(point[0],point[1],point[2])

    point_coords_data = doc.createTextNode(string)
    point_coords.appendChild(point_coords_data)

    #### Cells ####
    cells = doc.createElementNS("VTK", "Cells")
    piece.appendChild(cells)

    # Cell locations
    cell_connectivity = doc.createElementNS("VTK", "DataArray")
    cell_connectivity.setAttribute("type", "Int32")
    cell_connectivity.setAttribute("Name", "connectivity")
    cell_connectivity.setAttribute("format", "ascii")        
    cells.appendChild(cell_connectivity)

    # Cell connectivity data
    string = "\n"
    for i in range(intspace.space.simplices.shape[0]):
        string += " {:12d} {:12d} {:12d} {:12d}\n".format(
                    intspace.space.simplices[i,0],
                    intspace.space.simplices[i,1],
                    intspace.space.simplices[i,2],
                    intspace.space.simplices[i,3])
    connectivity = doc.createTextNode(string)
    cell_connectivity.appendChild(connectivity)

    cell_offsets = doc.createElementNS("VTK", "DataArray")
    cell_offsets.setAttribute("type", "Int32")
    cell_offsets.setAttribute("Name", "offsets")
    cell_offsets.setAttribute("format", "ascii")                
    cells.appendChild(cell_offsets)

    string1 = "\n"
    string2 = "\n"
    off    = 4
    for i in range(intspace.space.simplices.shape[0]):
        string1 += " {:12d}\n".format(off) # Offset
        string2 += " {:12d}\n".format(10)  # Type
        off += 4
    offsets = doc.createTextNode(string1)
    cell_offsets.appendChild(offsets)

    cell_types = doc.createElementNS("VTK", "DataArray")
    cell_types.setAttribute("type", "UInt8")
    cell_types.setAttribute("Name", "types")
    cell_types.setAttribute("format", "ascii")                
    cells.appendChild(cell_types)
    types = doc.createTextNode(string2)
    cell_types.appendChild(types)

    #### Data at Points ####
    point_data = doc.createElementNS("VTK", "PointData")
    piece.appendChild(point_data)

    point_data_id = doc.createElementNS("VTK", "DataArray")
    point_data_id.setAttribute("type", "Int32")
    point_data_id.setAttribute("Name", "PointID")
    point_data_id.setAttribute("format", "ascii")                
    point_data.appendChild(point_data_id)

    string = "\n"
    for i in range(intspace.space.points.shape[0]):
        string += " {:12d}\n".format(i)
    pointids = doc.createTextNode(string)
    point_data_id.appendChild(pointids)

    #### Data in Cells ####
    cell_data = doc.createElementNS("VTK", "CellData")
    piece.appendChild(cell_data)

    cell_data_id = doc.createElementNS("VTK", "DataArray")
    cell_data_id.setAttribute("type", "Int32")
    cell_data_id.setAttribute("Name", "SimplexID")
    cell_data_id.setAttribute("format", "ascii")                
    cell_data.appendChild(cell_data_id)

    string = "\n"
    for i in range(intspace.space.simplices.shape[0]):
        string += " {:12d}\n".format(i)
    cellids = doc.createTextNode(string)
    cell_data_id.appendChild(cellids)

    # Write to file and exit
    outFile = open(fname, 'w')
#        xml.dom.ext.PrettyPrint(doc, file)
    doc.writexml(outFile, newl='\n')
    outFile.close()
