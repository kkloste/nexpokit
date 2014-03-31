#!/bin/bash

# to call this shell script, use this:
# ./convert_data.sh
clear
echo Nexpokit dataset preparation

# change permissions to ensure this script can run everything it's supposed to
chmod 744 *.sh
chmod 744 *.m

echo dblp
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r convert_dblp > dblp.txt &

echo itdk
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r convert_itdk > itdk.txt &

echo flickr
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r convert_flickr > flickr.txt &

echo twitter
nohup /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r convert_twitter > twit.txt &

echo All conv files called!