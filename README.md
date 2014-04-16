FALCON
====

Flux Assignment (with) LAD Convex Objectives (and) Normalization.
LAD, or Least Absolute Deviations, is the 1-norm equivalent of least squares.



To INSTALL
====

Most users can directly build the minDisj program from the generated C file.

Users interested the minDisj file will need to install the ATS language, which
has a style similar to StandardML.

Once minDisj is built, consider adding the following to your startup.m
file (in UNIX, this should be in ~/Documents/MATLAB/), or run it
before each use of FALCON:

addpath(genpath('MY_FALCON_DIRECTORY'));

This will be especially useful if you plan on running analyses 
similar to those found in the publication; otherwise, the top-level
directory should suffice.


For other data and information related to FALCON and and associated
pubications, please see:

http://openwetware.org/wiki/Barker:Notebook/FALCON
