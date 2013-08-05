---
title: "Network Matrix Exponentials for link-prediction, centrality, and more"
layout: project
---


NEXPOKIT: Network Matrix Exponentials for link-prediction, centrality, and more
===============================================================================

### Kyle Kloster
### David F. Gleich

_These are research codes and may not work for you._

Download
--------

* [nexpokit.tar.gz](nexpokit.tar.gz) (2013-08-05)

Synopsis
--------

    compile % compile the mex files
    G = load_graph('dolphins');
    P = normout(G)';
    x = gexpmq_mex(P,1,11,1e-5,10*size(P,1));
    
Reusable codes
--------------

* `gexpm` a pure matlab prototype implementation of the Gauss-Southwell
  code with a heap
* `gexpm_mex` a C++ mex implementation of Gauss-Southwell with a heap
* `gexpmq_mex` a C++ mex implemetation of pseudo-Gauss-Southwell with
  a queue.
* `kmatexp` a matlab implementation of an N+1 step Taylor rule


Codes from others
-----------------

* `mexpv` from expokit
* `expv` from expokit

Results from the paper
----------------------

To reproduce figure 2 (left), run:

    test_tol_accuracy % generate the data
    plot_tol_accuracy % plot the data

To reproduce figure 2 (right), run:

    test_steps_accuracy_order
    plot_steps_accuracy
    
To reproduce figure 3, run

    test_runtime_tol4
    plot_runtime
    
    
        

Todo
----

1. Finish the readme
2. Implement the experiments
3. Mex the code

