% combine the data from the smaller graphs with webbase, twitter, and friendster
% /p/matlab-7.14/bin/matlab -nodisplay -nodesktop -nojvm -nosplash -r runtime_process > /scratch2/dgleich/kyle/joblog/runtime_proc.txt

experimentname = 'degruntime';
experiment_directory = '/scratch2/dgleich/kyle/nexpokit/results/';

degreedataname = 'degrees';

addpath('~/nexpokit');
%load(strcat('~/nexpokit/results/' , experimentname));
load(strcat(experiment_directory, experimentname));
load(strcat(experiment_directory, degreedatname)); % all the node-degree information is stored seperately (by accident)

% num_graphs
% num_degrees
% num_algs 

% time_vals = zeros(num_algs,num_degrees,num_graphs);
% err_vals = zeros(num_algs,num_degrees,num_graphs);
% alldegrees = zeros(num_degrees, num_graphs);

num_graphs = size(datasizes,1);

alldegrees = degrees;
errors = err_vals;
times = time_vals;
graphsizes = datasizes;

numrecords = num_graphs;

  % now load webbase
  	load(strcat(experiment_directory , experimentname, '_webbase'));
  	numrecords = numrecords + 1;
  	newdegreesnum = size(err_vals,2);
  	errors(:,1:newdegreesnum,numrecords) = err_vals(:,:);
  	times(:,1:newdegreesnum,numrecords) = time_vals(:,:);

	load(strcat(experiment_directory, degreedataname, '_webbase'));  %loads degrees  	
  	alldegrees(:,numrecords) = degrees(:);  	
  	graphsizes(numrecords,1) = datasizes(:);
  	
  % now load twitter
  	load(strcat(experiment_directory , experimentname, '_twitter'));
  	numrecords = numrecords + 1;
  	newdegreesnum = size(err_vals,2);  	
  	errors(:,1:newdegreesnum,numrecords) = err_vals(:,:);
  	times(:,1:newdegreesnum,numrecords) = time_vals(:,:);

	load(strcat(experiment_directory, degreedataname, '_twitter'));  %loads degrees  	
  	alldegrees(:,numrecords) = degrees(:);   	
  	graphsizes(numrecords,1) = datasizes(:);
  
  % now load friendster
  	load(strcat(experiment_directory , experimentname, '_friendster'));
  	numrecords = numrecords + 1;
  	newdegreesnum = size(err_vals,2);  	
  	errors(:,1:newdegreesnum,numrecords) = err_vals(:,:);
  	times(:,1:newdegreesnum,numrecords) = time_vals(:,:);

	load(strcat(experiment_directory, degreedataname, '_friendster'));  %loads degrees  	
  	alldegrees(:,numrecords) = degrees(:);   	
  	graphsizes(numrecords,1) = datasizes(:);	
	
% errors ( num_algs, num_degrees, num_data )
% times ( num_algs, num_degrees, num_data )
% graphsizes ( num_data, 1 )
% alldegrees ( num_degrees, num_data )


% HAVE ALL DATA for this experiment -- save in plotting/

save(strcat(experiment_directory, experimentname, '_to_plot', '.mat'), 'alldegrees', 'errors', 'times', 'graphsizes','-v7.3');
	
exit