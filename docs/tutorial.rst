graphy tutorial
==============

Example
-------------------------

Let's find a high-modularity decomposition of the karate-club network.

:doc:`graphy.partitions` implements search over decompositions.  After getting
a matrix, we need to create a cost function object (a subclass of from 
:doc:`graphy.costfunctions.CostFunction`), then pass it into 
:doc:`graphy.partitions.FindOptimal`.


.. plot:: test_pyplots/karate.py
   :include-source:


