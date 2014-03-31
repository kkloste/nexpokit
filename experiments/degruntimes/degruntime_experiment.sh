#!/bin/bash

# to call this shell script, use this:
# ./degruntime_experiment.sh
clear
echo Nexpokit degruntime experiment

# change permissions to ensure this script can run everything it's supposed to
chmod 744 *.sh
chmod 744 *.m

echo 'Begin small dataset trials output in   degrunsmall\n'
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r degruntime_experiment > degrunsmall.txt &

echo 'Begin twitter trials output in   degruntwit\n'
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r degruntime_experiment_twitter > degruntwit.txt &

echo 'Begin friendster trials output in   degrunfri\n'
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r degruntime_experiment_friend > degrunfri.txt &

echo 'Begin webbase trials output in   degrunweb\n'
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -singleCompThread -r degruntime_experiment_web > degrunweb.txt &


echo 'finished calling all experiments\n'