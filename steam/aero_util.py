""" Aerodynamic Utility Module

Module of common utilities for aero use in STEAM tool.
"""
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import steam
import pandas as pd
import numpy  as np
import math

# Degree trigonometry
def rad(x):
  return math.radians(x)

def deg(x):
  return math.degrees(x)

def cosd(x):
  return math.cos(rad(float(x)))

def sind(x):
  return math.sin(rad(float(x)))

def tand(x):
  return math.tan(rad(float(x)))

def acosd(x):
  return deg(math.acos(float(x)))

def asind(x):
  return deg(math.asin(float(x)))

def atand(x):
  return deg(math.atan(float(x)))

def atan2d(y,x):
  return deg(math.atan2(float(y),float(x)))

def atphi_to_ab(talpha,phi):
    """ Calculate alpha/beta from alphaT/phi.
    
    All values are in degrees.

    Returns (alpha,beta)
    """

    alpha = atan2d(cosd(phi)*sind(talpha),cosd(talpha))
    beta  = asind(sind(phi)*sind(talpha))

    return (alpha,beta)

def ab_to_atphi(alpha,beta):
    """ Calculate alphaT/phi from alpha/beta.
    
    All values are in degrees.

    Returns (alphaT,phi)
    """

    alphaT = acosd(cosd(alpha)*cosd(beta))
    phi    = atan2d(sind(beta),sind(alpha)*cosd(beta))

    return (alphaT, phi)


def fomoco(soln,A=1.0,L=1.0,mrc=None,weights=None):
    """ Integrate F/M on a solution that contains cp.

    Args:
        soln: A solution with a 'cp' variable
        A   : Reference Area         [1.0]
        L   : Reference Length       [1.0]
        mrc : Array for MRC location [0.0,0.0,0.0]
        weights: If you have the weights for this grid, then provie them (super fast)

    Returns:
        (:obj:`tuple`): array of forces, array of moments 
    """

    if (type(soln) != type(steam.solution.Solution())):
        raise TypeError("Input soln is not type(soln)!")

    if soln.store_at.upper().find('ELEM') != 0:
        raise Exception("Solution must be element-based!")

    if (weights is None):
        weights = fomoco_weights(soln.mesh)

    try:
        cp_vals   = soln.data.cp.values
    except:
        raise Exception("Could not get cp values!")

    fomoco = (cp_vals * weights).sum(axis=1) / A
    force  = fomoco[0:3]
    moment = fomoco[3:6]

    # Now do mrc converstion
    if mrc is not None:
        moment[0] = moment[0] + mrc[2] * force[1]  - mrc[1] * force[2]
        moment[1] = moment[1] - mrc[2] * force[0]  + mrc[0] * force[2]
        moment[2] = moment[2] + mrc[1] * force[0]  - mrc[0] * force[1]
    moment /= L

    return (force,moment)

@steam.util.timer(logger)
def fomoco_weights(mesh):
    """ Return a weight array for each element in the mesh for F/M integration.

    MRC is about [0,0,0]. Can be translated in fomoco call.

    Args:
        mesh: Mesh to integrate on

    Returns:
        (:obj:`array`): array of weights (nElem x 6)
    """

    if (type(mesh) != type(steam.mesh.Mesh())):
        raise TypeError("Input mesh is not type(mesh)!")

    mrc = np.zeros(3)

    #! Loop over every triangle and get the normal and area
    elements  = mesh.conn.values
    nodes     = mesh.xyz_pt.values
    elems     = mesh.xyz_el.values

    weights = np.zeros((6,len(elems)))

    for (e,elem) in enumerate(elements):
        n1 = nodes[elem[0]]
        n2 = nodes[elem[1]]
        n3 = nodes[elem[2]]
        nc = elems[e]
        #nc = (n1 + n2 + n3 )/ 3.0
        #print(n1,n2,n3,nc)

        v1 = n2 - n1
        v2 = n3 - n2
        v3 = np.cross(v1,v2)
        areas = v3 / 2.0
        #print(areas)

        fo = areas  * -1.0 # (Acts inward)
        mo = np.zeros(3)

        dx = nc - mrc

        mo[0] = fo[2] * dx[1] - fo[1] * dx[2]
        mo[1] = fo[0] * dx[2] - fo[2] * dx[0]
        mo[2] = fo[1] * dx[0] - fo[0] * dx[1]

        weights[0][e] = fo[0]
        weights[1][e] = fo[1]
        weights[2][e] = fo[2]
        weights[3][e] = mo[0]
        weights[4][e] = mo[1]
        weights[5][e] = mo[2]


    return weights


class LineLoads():
    """ Line loads object of saved weights/slices.

    This is provided to speed up processing of several
    solutions.  The weights and component centroids
    need only computed once and then can be rapidly applied to
    other solutions.

    Since line loads and possible debuggin involves a number
    of weights and a sliced mesh, a new class seemed prudent.

    This also allows creation of dedicated I/O for line loads.
    """
    
    def __init__(self,mesh=None,
                 A=1.0,
                 L=1.0,
                 res=1.0,
                 mrc=[0.0,0.0,0.0],
                 dir='X',
                 round=False,
                 ):
        """ Create LineLoads object.

        Args:
            mesh (:obj:`steam.mesh.Mesh()`): Mesh to cut, used by solutions
            res (:obj:`float`): resolution (in grid units) [1.0]
            A (:obj:`float`)   : Reference Area  [1.0]
            L (:obj:`float`)   : Reference Length  [1.0]
            mrc : Array for MRC location [0.0,0.0,0.0]
            dir : Direction of line loads ('X','Y','Z') ['X']
            round (:obj:`bool`): Round grid extents to integer numbers?  [F]

            A, L, and MRC can be overwritten in call to evaluate.

            Only 'X' direction is implimented (11/01/17)
        """

        if mesh is None:
            logger.error("LineLoads object needs a mesh defined")
            return

        self.mesh_orig  = mesh
        self.A          = A
        self.L          = L
        self.resolution = res
        self.mrc        = mrc
        self.dir        = dir.upper()

        #! Init the weights and arrays
        self.weights(self.mesh_orig,round)

    @steam.util.timer(logger)
    def weights(self,mesh,round=False):
        """ Returns 

        Args:
            mesh: Mesh to integrate on
            round (:obj:`bool`): Round to integer numbers?  [F]

        Returns:
            (:obj:`array`): array of weights (nElem x 6)
        """

        import math
        if (type(mesh) != type(steam.mesh.Mesh())):
            raise TypeError("Input mesh is not type(mesh)!")

        #! Operate on a copy of the mesh
        mesh = steam.mesh.Mesh(copy=self.mesh_orig)

        xmin = mesh.xyz_pt[self.dir].values.min()
        xmax = mesh.xyz_pt[self.dir].values.max()
        if round:
            xmin = math.floor(xmin)
            xmax = math.ceil( xmax)
        xrange = xmax - xmin

        ncomps  = math.ceil(
                            xrange / 
                            self.resolution
                            )  # Number of new components (buckets)
        nplanes = ncomps - 1                # Number of slice planes

        planes = []
        comps  = []
        cents  = []

        # Make a list of planes and new components.  The list of planes
        # are the cut locations.  On either side of the cuts components
        # are set to increment 1->ncomps
        for i in range(nplanes):
            x = xmin + (i+1) * self.resolution
            if self.dir == "X":
                planes.append([(x,0,0),(x,1,0),(x,0,1)])
            if self.dir == "Y":
                planes.append([(1,x,0),(0,x,0),(0,x,1)])
            if self.dir == "Z":
                planes.append([(1,0,x),(0,1,x),(0,0,x)])
            cents.append(xmin + (i+0.5) * self.resolution)
            if (i == 0):
                        comps.append([i+2,1])
            else:
                        comps.append([i+2,-1])
        # Add the last one
        cents.append(xmin + (nplanes+0.5) * self.resolution)


        # If nplanes = 0, then all of the grid fits inside of 
        #  one sliver of mesh.
        if (nplanes > 0):
            origElem = mesh.slice(planes,comps,returnList=True)
        else:
            origElem = np.array([i for i in range(mesh.conn.shape[0])])
            mesh.comp[:] = 1

        #! Ok... so I have the original element that it came from,
        #! so what I want is weightings from the original elements
        #! into the slice buckets.  The new components should give
        #! me the slice buckets.  The center of the buckes is known

        weights = steam.aero_util.fomoco_weights(mesh)
        #>>> weights.shape -> (6, NEL)

        #origElem ties me to the original mesh/soln.
        #weights are for the new elements' contribution to aero

        load_elems   = []
        load_weights = []

        for i in range(ncomps):
            icomp = i + 1  # components are 1-based
            ele = origElem[mesh.get_comp(icomp)]
            wei = weights[:,mesh.get_comp(icomp)]
            load_elems.append(ele)
            load_weights.append(wei)

        #! Store all of the outputs
        self.elements     = load_elems      # Original elements corresponding to weights
        self.weights      = load_weights    # Weights from each new element to the slice
        self.centroids    = np.array(cents) # Center of each slice

        # For debugging:
        self.mesh         = mesh       # Final mesh
        self.fomo_weights = weights    # In case we want to check later, save it now
        self.origElem     = origElem

        return

    #@steam.util.timer(logger)
    def calc(self,soln,A=None,L=None,mrc=None):
        """ Perform line load integration on provided solution.

        Arguments for A, L, and mrc default to what is in the
        object, but can be overwritten here.

        Args:
            soln               : Solution based on previous mesh
            A (:obj:`float`)   : Reference Area    [None]
            L (:obj:`float`)   : Reference Length  [None]
            mrc                : Array for MRC location [None]
        
        Returns:
            Numpy array of f/m coefficients on each slice
        """

        # Do the integration
        try:
                cp_vals   = soln.data['cp'].values
        except:
                raise Exception("Could not get cp values!")

        if mrc is None:
            mrc = self.mrc
        else:
            self.mrc = mrc

        if A is None:
            A   = self.A
        else:
            self.A = A

        if L is None:
            L   = self.L
        else:
            self.L = L
    
        fomocos = []
        for i in range(len(self.elements)):
            fomoco = (cp_vals[self.elements[i]] * self.weights[i]).sum(axis=1)
            fomocos.append(fomoco)
        fm = np.array(fomocos)

        fm /= A
 
        # Now do mrc converstion
        fm[:,3] = fm[:,3] + mrc[2] * fm[:,1]  - mrc[1] * fm[:,2]
        fm[:,4] = fm[:,4] - mrc[2] * fm[:,0]  + mrc[0] * fm[:,2]
        fm[:,5] = fm[:,5] + mrc[1] * fm[:,0]  - mrc[0] * fm[:,1]
        fm[:,3:6] /= L

        self.fomoco = fm # Save most recent
        return fm

    def write(self,name="output.ll",fm=None,header=[],fmt='%.12e'):
        """ Write the lineloads to a file.

        Args:
            name: Name of file to write [output.ll]
            fm  : optional output from LineLoad.calc() [None]
            header: optional lines to add to the header, array (no #)
        fmt (:obj:`string`): Format string according to numpy.savetxt(). ['%.12e']

        If fm not provided, then use the most recently calculated.

        """

        if fm is None:
            try:
                fm = self.fomoco
            except:
                raise Exception("Need to calculate line loads before you can write it!")

        #! Turn them into dataframes to dump with a cdat
        cent = pd.DataFrame(self.centroids,columns=[self.dir+"_Centroid"])
        fmdf = pd.DataFrame(fm,columns=["F_X","F_Y","F_Z","M_X","M_Y","M_Z"])
        cent = cent.join(fmdf)

        logger.info("Writing LL into cdat file '{}'".format(name))

        head  = ["Line load output file from STEAM",
                 "  The line load resolution is {} grid units.".format(self.resolution),
                 "  Reference Area: {:.9g} gu^2".format(self.A),
                 "  Reference Length: {:.9g} gu".format(self.L),
                 "  Moment Ref Center: [{:.9g}, {:.9g}, {:.9g}] gu".format(*self.mrc)]
        head.extend(header)
        vars  = cent.columns
        head.append("        ".join(vars))
        
        head_str = "\n".join(head)

        np.savetxt(name,cent.values,header=head_str,fmt=fmt)

    #@steam.util.timer(logger)
    def check(self,soln,A=None,L=None,mrc=None):
        """ Sanity checks for the line load.

        Compare result to the total integrated F/M as a sanity check.

        Args:
            soln               : Solution based on previous mesh
            A (:obj:`float`)   : Reference Area    [None]
            L (:obj:`float`)   : Reference Length  [None]
            mrc                : Array for MRC location [None]
        Returns:
            :obj:`bool`: Do things agree?  [T/F]

        """

        if mrc is None:
            mrc = self.mrc

        if A is None:
            A   = self.A

        if L is None:
            L   = self.L

        #! Does it compare to integrating the whole thing with the weights?
        fm = self.calc(soln,A,L,mrc)
        fm_ll = np.sum(fm,axis=0)
        
        soln2 = steam.solution.uniform_soln(self.mesh)
        soln2.add_var("cp",soln.data['cp'].values[self.origElem])
        # soln2.add_var("oID",origElem)  # <- Store old element IDs?
        fm_we = steam.aero_util.fomoco(soln2,self.A,self.L,self.mrc,weights=self.fomo_weights)

        # Check to make sure that they are all "very" close
        tol = 1e-10
        diff = abs(fm_ll[0:3] - fm_we[0])
        if (np.any(diff > tol)):
            logger.error("Line Loads forces  do not match entire integration")
            logger.error("{} {} {} > {}".format(*diff,tol))
            return False
        else:
            logger.info("Line Loads passed force check")

        diff = abs(fm_ll[3:6] - fm_we[1])
        if (np.any(diff > tol)):
            logger.error("Line Loads moments do not match entire integration")
            logger.error("{} {} {} > {}".format(*diff,tol))
            return False
        else:
            logger.info("Line Loads passed moment check")

        return True
