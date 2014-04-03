Nexpokit Experiments
===============================================================================

### Kyle Kloster
### David F. Gleich

** All parameters for all experiments are controlled entirely from [experiment name]_paramaters.m
For example, to set the number of trials used in the degruntime experiment to 57, simply open
'degruntime_parameters.m' and change "num_trials = 57". Other parameters: tol, t, maxnnz (for expm_svec)

Degree-based Runtime
-----------------

To execute all experiments, process all data, and plot:

(1) run './degruntime_experiment.sh'  (might require you to modify permissions)
(2) run './degruntime_plot.sh'
