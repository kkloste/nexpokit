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

To reproduce figures 1 and 2, run:
	
		experiment/localization_demo/example_localization_ljournal.m

To reproduce figure 4, first generate the data by running:

		experiments/accuracy_vs_error/compute_tol_accuracy_data.m

Then, to produce the plots for figure 4, run:

		experiments/accuracy_vs_error/plot_tol_accuracy.m
		
To generate the data for figure 5, run:

		experiments/accuracy_vs_work/compute_steps_accuracy_order.m

To reproduce figure 5, run:

		experiments/accuracy_vs_work/plot_tol_steps_accuracy.m
		
To reproduce figure 6, first generate the data by running:

		experiments/acc_vs_maxnnz/maxnnz_experiment.m
		experiments/acc_vs_maxnnz/maxnnz_experiment_web.m
		experiments/acc_vs_maxnnz/maxnnz_experiment_friend.m				
		experiments/acc_vs_maxnnz/maxnnz_experiment_twitter.m

Then, to produce the plots for figure 6, run

		experiments/acc_vs_maxnnz/maxnnz_plots_finalized.m

To reproduce figure 7, first generate the data by running:

		experiments/runtimes/runtime_experiment.m
		experiments/runtimes/runtime_experiment_web.m		
		experiments/runtimes/runtime_experiment_friend.m		
		experiments/runtimes/runtime_experiment_twitter.m
		experiments/runtimes/runtime_process.m		

Then, to produce the plots for figure 7, run:

		experiments/runtimes/runtime_plot.m

To reproduce figure 8, first generate the data by running:

		experiments/scaling/scaling_study_1.m
		experiments/scaling/scaling_study_s.m

Then, to produce the plots for figure 8, run:

		experiments/scaling/scaling_plots.m

