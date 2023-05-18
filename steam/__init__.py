""" Autoload modules in the STEAM package. """

__version__ = "2.0"

#! Logging Initialization
#! - Do it before we load all the other imports so that they can baseline this
#!     https://docs.python.org/3/howto/logging.html#logging-basic-tutorial

import logging, logging.config

logLevel_file  = 'DEBUG'
logLevel_screen = 'WARN'
logConfig = {
            'version':1,
            'disable_existing_loggers':False,
            'formatters':{
                'default':{
                    'format':'%(asctime)s  %(levelname)-5s  %(name)-20s: %(message)s',
                    },
                'screen':{
                    'format':'%(levelname)-5s  %(name)-20s: %(message)s',
                    },
            },
            'handlers':{
                # This is the logger for STEAM
                'steam_file':{
                    #'class':'logging.handlers.RotatingFileHandler',
                    #'maxBytes':10*1024*1024,  # 10M size for logs
                    #'backupCount':0,
                    'class':'logging.FileHandler',
                    'formatter':'default',
                    'filename':'steam.log',
                    'mode':'w',            # Replace with new file
                    'level':logLevel_file,
                },
                'steam_screen':{
                    'class':'logging.StreamHandler',
                    'formatter':'screen',
                    'level':logLevel_screen,
                },
            },
            'loggers':{
                # This is the root logger and captures everything
                '':{
                    'handlers':['steam_file','steam_screen'],
                    'propagate':True,
                    'level':logLevel_file,
                },
            },
         }
logging.config.dictConfig(logConfig)

logger = logging.getLogger(__name__)

def log_level(file=None,screen=None):
    """ Change the file log level."""

    for hand in logging.getLogger().handlers:
        if hand._name == "steam_file":
            log_hand_file   = hand
        if hand._name == "steam_screen":
            log_hand_screen = hand

    if file  is not None and log_hand_file is not None:
        logger.debug("Changing file log level to {}".format(file))
        log_hand_file.setLevel(file)

    if screen is not None and log_hand_screen is not None:
        logger.debug("Changing screen log level to {}".format(screen))
        log_hand_screen.setLevel(screen)

def log_off(file=True,screen=True):
    """ If true, will deactivate that log.
    
        Defaults to deactivating all logging.
    """
    if screen:
        logger.debug("Deactivating screen log")
        change_log_level(level_screen=100)
    if file:
        logger.debug("Deactivating file log")
        change_log_level(level_file=100)
    
# Done with logging set-up


logger.info("Importing Subpackages")

# Import all other modules
from . import mesh
from . import solution
from . import table
from . import database
from . import interpolate
from . import io
from . import util
from . import aero_util
from . import container
from . import models
logger.info("Done importing Subpackages")

#! Only import SWIG items if they are compiled
#! init.py is a proxy test for all of them.
has_libmesh = False
import os
dir_root = os.path.dirname(os.path.abspath(__file__))
if os.path.isfile(dir_root+'/libmesh/init.py'):
    logger.info("Importing LibMesh")
    has_libmesh = True
    from . import libmesh
    from . import physics_models
    logger.info("Done importing LibMesh")

#! Try to import the PyData package (import pydata)
#! This is a separate package which can be installed
#! elsewhere and is accessed by adding the path to
#! the package to the users PYTHONPATH.
#! export PYTHONPATH=/path/to/pydata/:$PYTHONPATH
has_dataio = False
try:
    import pydata
    has_dataio = True
    logger.info("Done importing PyData")
except ImportError:
    pass

