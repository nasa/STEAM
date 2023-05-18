Mesh
----------

CEbner: Conceptually, what is a mesh?  Describe the two different types (point cloud versus ?not? point cloud). What
is our parlance for the three: node-based mesh, element-based mesh, and point cloud?

A mesh is the fundamental unit in which geometry is stored in STEAM.  STEAM
meshes are surface only; no volume meshes are used and therefore volume data
is not used.  Meshes are fundamentally three dimensional and unstructured.
Mesh elements, therefore, are triangular with three nodes associated and 
up to three adjacent elements.  (Elements can have fewer than three neighbors
if they reside along an edge of the mesh.)

CEbner: How is a mesh generally created?  Is it made whole cloth or imported from a file?  What filetypes are
supported or where can I look to find them?

CEbner: What are the user-facing attributes of interest and what is their form?  Is it different if the mesh is
a point cloud or a non-point cloud.  How do the attribures change if we decide to make it  node-based or element-based?

CEbner: What is a static solution and how do I create one?  How do I write it out?


Components
~~~~~~~~~~~~~~~~~~~~

CEbner: What are components (index-based and names)?

CEbner: What are some common reasons to have components?  Do they have impact on solutions?

CEbner: How do I create components for my mesh?

Simple Creation
~~~~~~~~~~~~~~~~~~~~
CEbner: How can I create a mesh to play around?  STEAM can make spheres and cubes easilly.  Show an example?

.. _mesh_point_cloud
Point Cloud
~~~~~~~~~~~~~~~~~~~~

A Point Cloud is a specific type of Mesh with no connectivity.  This is useful
for capturing, storing, and interacting with data in three dimensional space
but without elements.

CEbner: Why would I want to, or have to, use a point cloud instead of a mesh?  Are those the correct terms?

Examples
~~~~~~~~~~~~~~~~~~~~
Files referenced in this section are collected in the `examples/` directory.  To not clutter 
the installation directory, copy the `examples/` directory and use it as your working directory 
for the following examples.

.. contents:: :local:

Read / Write
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: Example of reading and writing a mesh.  Perhaps create a mesh in STEAM, write to format A, read from file, write to format B.


Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: We can perform transformations, is that done in-place?  Can you show an example?

Node to Element
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: How can I perform a node-to-element conversion?

Components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CEbner: Can we provide an example of reading components (or using the ones on mk_sphere) to create nested components?
CEbner: Can we extract specified components into a new mesh? 
CEbner: Show that the solution information is transfered with it.
