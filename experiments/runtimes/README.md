Nexpokit Experiments
===============================================================================

### Kyle Kloster
### David F. Gleich

** All parameters for all experiments are controlled entirely from [experiment name]_experiment.m
For example, to set the number of trials used in the runtime experiment to 57, simply open
'runtime_experiment.m' and change "num_trials = 57". Other parameters: tol, t, maxnnz (for expm_svec)

Runtime
-----------------

To execute all experiments, process all data, and plot:

(1) run 'runtime_experiment.sh' using nohup ./runtime_experiment.sh > expstatus.txt &, and once it's done, then
(2) run 'runtime_plot.sh' using ./runtime_makeplot.sh > plotstatus.txt &

Top-k accuracy
----------------------
