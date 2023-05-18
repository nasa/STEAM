.. DBMS documentation master file, created by
   sphinx-quickstart on Wed Jan 25 09:51:10 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the STEAM Documentation
==================================

For an introduction into usage, please refer to the sections of the User's Guide:

.. toctree::
   :maxdepth: 2

   users_guide

.. toctree::
   :maxdepth: 1

   developers_guide


For developers and users interested in the STEAM API, see the following documentation:

.. toctree::
   :maxdepth: 2
   :caption: Python-based Modules:
   :glob:

   steam

   steam.aero_util
   steam.container
   steam.database
   steam.interpolate
   steam.io
   steam.mesh
   steam.solution
   steam.table
   steam.models
   steam.util

.. only:: libmesh

   .. toctree::
      :maxdepth: 2
      :caption: LibMesh-based Modules:
   
      steam.libmesh
      steam.physics_models


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
