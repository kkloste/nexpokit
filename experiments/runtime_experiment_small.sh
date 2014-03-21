#!/bin/bash

# to call this shell script, use this:
# ./runtime_experiment_small.sh
clear
echo Nexpokit runtime small experiment

# change permissions to ensure this script can run everything it's supposed to
chmod 744 *.sh
chmod 744 *.m

echo Begin small dataset trials output in   runsmall
/p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r runtime_experiment_small > runsmsmall.txt

echo finished calling all experiments

echo Collect the runtime data from all experiments and process it
echo Then make the plots
/p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nosplash -r runtime_plot_small > runsmsmall.txt

echo Plot made!