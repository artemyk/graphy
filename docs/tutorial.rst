graphy tutorial
==============

Example
-------------------------

Here we find a high-modularity decomposition of the karate-club network.

:doc:`graphy.partitions` implements search over decompositions.  After getting
a matrix, we need to create a cost function object (a subclass of from 
:doc:`graphy.qualityfuncs.QualityFunction`), then use the `find_optimal` method.
This returns a membership vector, 
which, for an N-dimensional system, is an N-dimensional integer-valued 
numpy array indicating the community each component belongs to.


.. plot:: test_pyplots/karate.py
   :include-source:



