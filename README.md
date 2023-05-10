# Overview

Welcome to the main code repository for the  **S**oftware 
**T**oolkit for **E**ngineering and **A**eroscience **M**odels (STEAM). STEAM
development began at NASA's Johnson Space Center (JSC) within the Applied 
Aerosciences and Computational Fluid Dynamics Branch (EG3)
with the goal of providing a toolbox enabling interpolation and 
model application in aerodynamic and aerothermodynamic databases of 
Computational Fluid Dynamics (CFD) solutions across varying geometry and flight
spaces. 

Using STEAM, one can assemble databases comprising CFD solutions from
multiple solvers on differing grids.  Geometric manipulation and interpolation 
tools are included to assimilate solutions onto common meshes and perform
simple mesh modification.
Flight-space interpolation can be done
using user-defined independent variables and subsequently applied models 
need not share the same independent parameters.  This flexibility allows a 
user to define the most useful databases and models for the problem at hand,
rather than being constrained by the software available to them.

### Operation in JSC's Flight Sciences Laboratory
STEAM was initially, and still primarily is, developed in the computing
environment of 
JSC's Flight Sciences Laboratory (FSL).  Included in this repository, 
therefore, are some files, scripts, and bits of documentation that may not
be valuable to others operating outside of the FSL.  Despite STEAM's 
inextricable link to the FSL, it runs on systems
outside of the FSL and is used in environments ranging from personal laptops
to large High Performance Computing (HPC) clusters.

# Git Branching Guidelines
Since its inception, STEAM has operated with a 
git branching strategy based on 
[this blog post](http://nvie.com/posts/a-successful-git-branching-model).
Please read at least the sections on 
[the main branches](https://nvie.com/posts/a-successful-git-branching-model/#the-main-branches) and 
[feature branches](https://nvie.com/posts/a-successful-git-branching-model/#feature-branches) and be sure to adhere to the guidance provided.

The main branches are **master** and **develop**.  Feature branches are 
traditionally prefixed by your last name in a 
directory-like structure.  An example of several feature branches might be:

* schwing/feature_1
* schwing/bugfix_1
* schwing/feature_2

Remember to merge features back to **develop**.  Merges back 
to **master** should be done for releases.  Further, almost no one should be
merging into master _ever_.  If you're not certain whether or not you're one of
the ordained few who can merge into master, you're not.
Any code that is checked into **develop** should pass the regression/V&V suite.

# STEAM's Python Environment
STEAM currently runs in an environment leveraging Python 3.10.  The development
team has traditionally maintained the necessary
Python environment using 
[Anaconda](https://docs.conda.io/en/latest/miniconda.html), though it isn't
strictly necessary for one with experience managing Python environments.
For those using Anaconda, an environment `py310.yml` file is provided in the 
repository to expedite the process of creating the proper environment.

Users in the FSL can further leverage the FSL's module system with the provided
STEAM module for easily establishing the proper versions of Anaconda, MPI, and
other executables.
