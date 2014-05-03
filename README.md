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

* `gexpm_mex` a C++ mex implementation of Gauss-Southwell with a heap
* `gexpmq_mex` a C++ mex implemetation of pseudo-Gauss-Southwell (rounded Gauss-Seidel) with
  a queue.
* `expmimv_mex` a C++ mex implemetation of an N step Taylor polynomial via Horner rule and "incomplete" matrix-vector products
  a queue.
* `kmatexp` a matlab implementation of an N+1 step Taylor polynomial via Horner rule
* `taydeg.hpp` : contains function for selecting Taylor degree appropriate for input error of “tol”
* `heap.hpp` : contains heap functions used in gexpm_mex


Codes from others
-----------------

* `mexpv` from expokit
* `expv` from expokit
For Higham & Al-Mohy's "expmv":
* `expmv`
* `expmv_tspan`
* `normAm`
* `select_taylor_degree
* `theta_taylor_half.mat`
* `theta_taylor_single.mat`
* `theta_taylor.mat`


Results from the paper
----------------------

To reproduce figure 2 (left), run:
	[fig 2 shows the non zeros of the solution, compared to the non zeros used by each algorithm for a particular input tolerance]
	
	experiment/localization_demo/example_localization_ljournal.m

[ Fig 4 is log10 of error vs. precision ]
To generate the data for figure 4, run:

		experiments/accuracy_vs_error/compute_tol_accuracy_data.m

To reproduce figure 4, run:

		experiments/accuracy_vs_error/plot_tol_accuracy.m
		
[ Fig 5 is work vs. top-k precision ]
To generate the data for figure 5, run:

		experiments/accuracy_vs_work/compute_steps_accuracy_order.m

To reproduce figure 5, run:

		experiments/accuracy_vs_work/plot_tol_steps_accuracy.m
		

[ Fig 6 is maxxnnz vs. error/precision ]
To reproduce figure 6, run:

		experiments/acc_vs_maxnnz/maxnnz_plots_finalized.m

[ Fig 7 is size vs. runtime ]
To reproduce figure 7, run:
runtime_plot.m

[ Fig 8 is synthetic experiments ]
scaling study 1
scaling study 1s
(generate raw data ^)
scaling plot (actually generates plot)
To reproduce figure 8 (left), run:
To reproduce figure 8 (middle), run:
To reproduce figure 8 (right), run:
