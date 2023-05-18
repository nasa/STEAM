Solution
--------
A solution is a repository for data mapped to a Mesh_.  Therefore, each solution
has an associated mesh, stored in the mesh attribute.
Similar to CFD applications, an arbitrary number of solutions can be 
associated with one mesh.

Solutions can hold data at either mesh nodes or mesh element centroids.

Solutions have several attributes which are almost always populated in use:

* **mesh:** The mesh with which the solution is associated.
* **store_at:** Is data stored at mesh nodes or mesh element centroids?  Possible values are "NODES" or "ELEMENTS".
* **data:** a pandas.DataFrame storing all of the mapped data.  The number of columns in **data** will equal the number of variables in the solution and the number of rows will equal either the number of nodes or elements in the mesh, depending on **store_at**.

CEbner: Based on what was filled-in for the Mesh section, does the above need updated?

CEbner: What is the difference between node- and element-based solutions?  How do I go between the two?

CEbner: Discuss how going node-to-element or back diffuses data and can reduce accuracy.

Simple Creation
~~~~~~~~~~~~~~~~~~~~

CEbner: What methods exist to create a dummy solution for examples and experiments?

Models (?)
~~~~~~~~~~~~~~~~~~~~

CEbner: Does this section need to exist?  Should it be here or somewhere else?

CEbner: What is a model?  What types of models are there and why do we even need to talk about them?

Examples
~~~~~~~~~~~~~~~~~~~~
Files referenced in this section are collected in the `examples/` directory.  To not clutter 
the installation directory, copy the `examples/` directory and use it as your working directory 
for the following examples.

.. contents:: :local:

Read / Write
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: How do I read and write a solution?  Format conversion?

Node to Element
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: Example of going from node to elemement or vise-versa.  Can we show an image on a cube corner?

Components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: How can I use components (of a mesh) to select a subset of a solution?


Models (?)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: Does this section need to exist?  Should it be here or somewhere else?

CEbner: Example of using the availible models that we have to augment solutions.

