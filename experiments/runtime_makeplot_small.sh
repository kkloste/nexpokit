#!/bin/bash

# to call this shell script, use this:
# ./runtime_makeplot_small.sh
clear
echo Nexpokit runtime small plot

# change permissions to ensure this script can run everything it's supposed to
chmod 744 *.sh
chmod 744 *.m

echo Collect the runtime data from all experiments and process it
echo Then make the plots
/p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r runtime_plot_small > runtime_plot.txt

echo Plot made!