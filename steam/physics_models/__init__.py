## We do a check here to make sure that SWIG is built
## Otherwise, it totally kills the Sphinx build
import logging
logger = logging.getLogger(__name__)
logger.debug("Top of file")

import os
dir_root = os.path.dirname(os.path.abspath(__file__))
if os.path.isfile(dir_root+'/modified_newtonian.py'):
    from . import modified_newtonian
    from . import locate_stagnation_point
    from . import streamlines

class StreamlineContainer:

    def __init__ (self,mesh):
        self.parent = mesh
        self.p = streamlines.get_new_pointer(mesh())

    def __del__ (self):
        streamlines.delete_pointer(self.p)

    def __call__(self):
        return self.p

    def create (self,stag_method,stream_method,fs):
        streamlines.create(self.p,stream_method,fs)

    def clear (self):
        streamlines.clear(self.p)

    def write (self,filename):
        streamlines.write(self.p,filename)
